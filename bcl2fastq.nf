#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         NF-bcl2fastq
========================================================================================
 NF-bcl2fastq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/gnetsanet/NF-bcl2fastq
 #### Authors
 Netsanet Gebremedhin gnetsanet <gnetsanet@gmail.com> - https://github.com/gnetsanet>
----------------------------------------------------------------------------------------
*/

import groovy.util.*
import java.io.File

/*
	Draft POC Nextflow script for Bcl2FASTQ processing
	Usage example: nextflow run /mnt/home/ngebremedhin/novaseq_pipe/hts-scripts/nextflow-novaseq/bcl2fastq.nf \
	--samplesheet /mnt/galaxy/novaseq_staging/171106_A00336_0009_AH57HFDMXX/SampleSheet.csv \
	-process.executor sge -with-dag -with-report -with-timeline -resume -N ngebremedhin@neb.com \
	--link false --nomail true --dryrun false

*/

// To-do:
// 1. DONE - sort of - Parse the experiment name from the samplesheet using NF or Groovy CSV parser. Not terribly important but is used in the email notification
// 2. Java memory allocation is hard-coded right now à la -Xmx64000000k - Not super urgent since we can specify resources per process.
// 3. DONE - still lots of assumptions - The number of rows to skip - in the SampleSheet.csv - is hard-coded. Should be changed to something dynamic. Particularly so if it the number of header rows varies

// Defaults and input params ****************************************************************************************************
samplesheet = params.samplesheet
params.max_mismatches = params['cutoffs_and_thresholds'].max_mismatches
params.max_no_calls = params['cutoffs_and_thresholds'].max_no_calls
params.min_mismatch_delta = params['cutoffs_and_thresholds'].min_mismatch_delta
params.extract_bc_cpus = params['cutoffs_and_thresholds'].extract_bc_cpus
params.tile_process_cpus = params['cutoffs_and_thresholds'].tile_process_cpus
params.picard = params['script_paths'].picard
params.java = params['script_paths'].java
params.destination = params['script_paths'].destination
params.post_run_scripts = params['script_paths'].post_run_scripts
params.fasqtc = params['script_paths'].fastqc
params.multiqc = params['script_paths'].multiqc
params.link = true // whether or not to create symlinks to combined fastqs. Set to false while testing since seq-shepherd, too does this.
params.nomail = false // Good parameter to have for testing since we do not want to spam people
params.mail_function = 'ruby'
params.dryrun = true
params.testmode = false // Perform end-to-end test with a test dataset - take a fraction of the data ( two tiles) from any instrument run
params.kraken = params['script_paths'].kraken
params.krakendb = params['script_paths'].krakendb
params.ktImportTaxonomy = params['script_paths'].ktImportTaxonomy
paras.run_kraken_classify = false // Do not run kraken classify by default - has to be explicity turned on.


//=============================================== Staging and Preprocessing ====================================================

metrics_name="${params.max_no_calls}nc_${params.min_mismatch_delta}mmd_${params.max_mismatches}mis_bc_metrics.txt"
runDir = getRunDirectory(samplesheet)
barcodesOutputDir = createOutputDir(samplesheet)
(tile_xml, instrument) = getTileXmlFile(samplesheet)
(tilesList, preProcessDir) = retrieveListOfTiles(tile_xml, instrument, runDir)
(num_reads,run_id, flowcell, machine_name, lanecount, read_structure) = parseRunXml(runDir)
(headerLines, experimentName) = parseSampleSheet(samplesheet)


// I do not think the generateBarcodeParams process is needed since the following block of code is meant to replace that
// And the following for loop ( over lanes) along with the samplesheet accomplishes what generateBarcodeParams is doing.
for( int laneID in 1 .. lanecount.toInteger()) {

	String laneNo = laneID
	File barcodeParamsFile = new File("${preProcessDir}/lane${laneNo}_barcode_params.txt")
	File multiplexParamsFile = new File("${preProcessDir}/lane${laneNo}_multiplex_params.txt")
	if (num_reads > 3) {
		barcodeParamsFile << "barcode_sequence_1\tbarcode_sequence_2\tbarcode_name\tlibrary_name\n"
		multiplexParamsFile << "OUTPUT_PREFIX\tBARCODE_1\tBARCODE_2\n"
		multiplexParamsFile << "L${laneNo}_unassigned\tN\tN\n"
		Channel.fromPath(samplesheet)
			.splitCsv(header: true, skip: headerLines, strip: true)
			.map { row ->
					barcodeParamsFile << "${row.index}\t${row.index2}\t${row.I7_Index_ID}+${row.I5_Index_ID}\t${row.Sample_ID}\n"
					multiplexParamsFile << "L${laneNo}_${row.Sample_ID}\t${row.index}\t${row.index2}\n"
			}
	} else {
			barcodeParamsFile << "barcode_sequence_1\tbarcode_name\tlibrary_name\n"
			multiplexParamsFile << "OUTPUT_PREFIX\tBARCODE_1\n"
			multiplexParamsFile << "L${laneNo}_unassigned\tN\n"
			Channel.fromPath(samplesheet)
				.splitCsv(header: true, skip: headerLines, strip: true)
				.map { row ->
						barcodeParamsFile << "${row.index}\t${row.I7_Index_ID}\t${row.Sample_ID}\n"
						multiplexParamsFile << "L${laneNo}_${row.Sample_ID}\t${row.index}\n"
				}
	}
}

// Convert text files to channels so they can be used as inputs in downstream processes
barcodeFiles = Channel.fromPath("${preProcessDir}/lane*_barcode_params.txt")

// Create an input that can be consumed by 'copyCombineFastq' - a combination of sample and read { 1, 2} to combined per tile fastqs into

Channel.fromPath(samplesheet)
  .splitCsv(header: true, skip: headerLines, strip: true)
	.map { row ->
		libraryName = row.Sample_ID.replaceAll("\\s","-") //Replace one or more spaces within library names
	}
	.unique()
	.into { sampleNames; sampleNamesForKraken; sampleNamesForKrona }

sampleNames
	.spread(Channel.from(getReadsToTransfer(num_reads)))
	.set { reads2Transfer }

lanes = Channel.from(1..lanecount)

// Create a list of  available lane and tile combinations that can be used to parallelize the bcl2fastq process
// tilesList is returned by a function that parses the run xml file
Channel.from(tilesList)
  .splitCsv(header: false, limit: params.testmode ? 2 : 0) // Counter-intuitively, limit of 0 means everything/No limit.
  .map {
  	lane= (instrument == 'novaseq') ? it[0].split("_")[0] : 1
  	tile= (instrument == 'novaseq') ? it[0].split("_")[1] : it[0]
  	[lane, tile]
  }
  .set { tiles }

// Two channel subsetting and transformation to generate inputs for use by linkFastqsToUsersFolders
// and notifyUsers processes

Channel.fromPath(samplesheet)
        .splitCsv(header: true, skip: headerLines, strip: true)
        .map { row ->
                sampleId = row.Sample_ID.replaceAll("\\s","-") //Strip a csv field's leading and trailing white spaces with splitCsv's 'strip : true' option
                def m = row =~ /[a-zA-Z0-9_.]+\@neb.com/; // project owner's email may not be in a specific column albeit it usually is in a column called 'Description'
								user = m[0]
                [sampleId, user]
        }
        .into { sampleAndUserFastqc; sampleAndUserLinkFastq }

Channel.fromPath(samplesheet)
  .splitCsv(header: true, skip: headerLines, strip: true)
	.map { row ->
		def m = row =~ /[a-zA-Z0-9_.]+\@neb.com/;
		m[0]
	}
	.unique()
	.into { user_emails_nextflow; user_emails_ruby }

//End of Staging and Preprocessing =======================================================================================================

// Start of Processes ####################################################################################################################

process ExtractIlluminaBarcodes {
	tag { ["ExtractIlluminaBarcodes", lane ] }
	publishDir "${barcodesOutputDir}", mode: 'move', overwrite: true

	when: !params.dryrun

	input:
	file barcodeFile from barcodeFiles
	val lane from lanes

	output:
	file("*_barcode.txt") into barcodesExtracted

	script:
	"""
	${params.java}/java -Xmx64000000k -jar ${params.picard}/picard.jar \
  ExtractIlluminaBarcodes \
	MAX_NO_CALLS=${params.max_no_calls} MIN_MISMATCH_DELTA=${params.min_mismatch_delta} \
	MAX_MISMATCHES=${params.max_mismatches} NUM_PROCESSORS=${params.extract_bc_cpus} \
	read_structure=${read_structure} \
	LANE=${lane} \
	BASECALLS_DIR="${runDir}/Data/Intensities/BaseCalls/" \
	METRICS_FILE="${barcodesOutputDir}/L${lane}_${metrics_name}" \
	BARCODE_FILE=${barcodeFile} \
	OUTPUT_DIR=.
	"""
}

barcodesExtracted.into { beginFlagBarcodeSummary; beginFlagBcl2Fastq }

process createBarcodeSummary {
	when: !params.dryrun

	input:
	file extractionComplete from beginFlagBarcodeSummary.toList()

	script:
	Channel.fromPath("${barcodesOutputDir}/L*_${metrics_name}")
		.first()
		.collectFile(skip: 6) // ignore header which is the first 6 lines
		.splitText()
		.filter({it.tokenize('\t').size() >= 3}) // split tab-delimited file by tab and ignore blank lines using the count of column data as a proxy for blankness.
		.map { row ->
			columns=row.tokenize("\t")

			if (columns[0] =~ /^N+$/) {
				columns[2] = 'Unknown' // rename columns since no library or sample name for unassigned barcodes and none is output in the file we are processing
				columns[3] = 'Undetermined'
			}
			[columns[0], columns[1], columns[2], columns[3], columns[4], columns[11]].join("\t")

		}
		.collectFile(
			name: "${barcodesOutputDir}/barcode_summary.txt", newLine: true, sort: false
		)
		'''
		'''
}

process bcl2fastq {
	tag { [lane, tile] }
	publishDir "${barcodesOutputDir}/fastq/L_${lane}_${tile}", mode: 'move', overwrite: true

	when: !params.dryrun

	input:
	set lane, tile from tiles
	file barcodesExtractCompletionMarker from beginFlagBcl2Fastq.toList()

	output:
	file("*.fastq.gz") into fastqFiles


	script:
	"""
  ${params.java}/java -Xmx16000000k -jar ${params.picard}/picard.jar IlluminaBasecallsToFastq \
  NUM_PROCESSORS=${params.tile_process_cpus} \
  read_structure=$read_structure \
  RUN_BARCODE=${machine_name} \
  LANE=${lane} \
  FIRST_TILE=${tile} \
  TILE_LIMIT=1 \
  MACHINE_NAME=${machine_name} \
  FLOWCELL_BARCODE=${flowcell} \
  BASECALLS_DIR="${runDir}/Data/Intensities/BaseCalls/" \
  BARCODES_DIR=${barcodesOutputDir} \
  MULTIPLEX_PARAMS="${preProcessDir}/lane${lane}_multiplex_params.txt" \
  MAX_READS_IN_RAM_PER_TILE=3000000 \
  MAX_RECORDS_IN_RAM=3000000 \
  COMPRESS_OUTPUTS=true \
  COMPRESSION_LEVEL=1 \
  TMP_DIR=/state/partition1/sge_tmp
	"""
}


process copyCombineFastq {

	tag { sample }
	publishDir "${barcodesOutputDir}" , mode: 'move', overwrite: true

	when: !params.dryrun

	input:
	set sample, read from reads2Transfer
	file fastqs from fastqFiles.toList()

	output:
	file("${sample}_${read}_Combine_Complete.txt") into combineFastqCompletionMarkers

	script:
	"""
	pushd ${barcodesOutputDir}
	sampleId=\$(echo ${sample} | sed -r 's/[: /\\]/-/g' | sed  's/^L._//' | tr "\\n" " " )
	/bin/bash ${params.post_run_scripts}/copy_combine_fastq.sh \${sampleId} ${read}
	popd
	touch ${sample}_${read}_Combine_Complete.txt
	"""
}

combineFastqCompletionMarkers.into { combineFastqCompletionMarkersForKraken; combineFastqCompletionMarkersForFastqc }


process runKrakenClassify {

	tag { sampleId }
  publishDir "${barcodesOutputDir}/kraken/" , mode: 'move', overwrite: true

	when: !params.dryrun && paras.run_kraken_classify

	input:
	val sampleId from sampleNamesForKraken
	file combinecompletionMarker from combineFastqCompletionMarkersForKraken.toList()

	output:
	file("${sampleId}.kraken") into krakenOut


	script:
		if(getReadsToTransfer(num_reads).size()>1) {
			"""
			${params.kraken} --db  ${params.krakendb} --paired ${barcodesOutputDir}/combined_fastq_tmp/${sampleId}.1.fastq.gz ${barcodesOutputDir}/combined_fastq_tmp/${sampleId}.2.fastq.gz --threads 4 > ${sampleId}.kraken
			"""
		} else {
			// It is better to create a new matrix spread that contains combination of library id and fastq paths than repeat ${barcodesOutputDir}/combined_fastq_tmp/
			"""
			${params.kraken} --db ${params.krakendb} ${barcodesOutputDir}/combined_fastq_tmp/${sampleId}.*.fastq.gz --threads 4 > r${sampleId}.kraken
			"""
		}

}


process createKronaInput {

	tag { library }
  publishDir "${barcodesOutputDir}/kraken/" , mode: 'move', overwrite: true

	when: !params.dryrun && paras.run_kraken_classify

	input:
	val library from sampleNamesForKrona
	file ("${library}.kraken") from krakenOut

	output:
	file("${library}.krona.in") into kronaIn

	script:
	"""
	cut -f2,3 ${library}.kraken > ${library}.krona.in
	"""
}


process runFastQc {

  tag { sampleId }
  publishDir "${barcodesOutputDir}/fastqc/fastqc_${sampleId}_logs" , mode: 'copy', overwrite: true

	when: !params.dryrun

	input:
  set sampleId, user from sampleAndUserFastqc // May not need info about the user here. Revisit: List of samples should suffice.
	file fastqs from combineFastqCompletionMarkersForFastqc.toList()

  output:
	// file("${sampleId}_FastQC_Complete.txt") into fastQcCompletionMarkers
	file("*_fastqc.zip") into fastQcCompletionMarkers

  script:
  """
  #mkdir -p ${barcodesOutputDir}/fastqc/fastqc_${sampleId}_logs
  ${params.fasqtc} -o ${barcodesOutputDir}/fastqc/fastqc_${sampleId}_logs -f fastq -q ${barcodesOutputDir}/combined_fastq_tmp/${sampleId}.*.fastq.gz
	#touch ${sampleId}_FastQC_Complete.txt
	"""
}


process generateKronaReport {

	tag { "Generating Krona html Report..."}
	publishDir "${barcodesOutputDir}/kraken/", mode: 'move', overwrite: true

	when: !params.dryrun && paras.run_kraken_classify

	input:
	file kronaInputs from kronaIn.toList()

	output:
	file('krona_report.html') into kronaReport


	script:
	"""
	${params.ktImportTaxonomy} ${kronaInputs} -o krona_report.html
	"""

}

process runMultiQc {
  tag { "multiQC"}
  publishDir "${barcodesOutputDir}", mode: 'move', overwrite: true

	when: !params.dryrun

	input:
	file fastqcMarker from  fastQcCompletionMarkers.toList()

  output:
	file('multiQC_Complete.txt') into multiQcCompletionMarker

  script:
  """
  pushd ${barcodesOutputDir}/fastqc
  ${params.multiqc} .  --force
  popd
	touch 'multiQC_Complete.txt'
  """

}


process linkFastqsToUsersFolders {
	tag {[user, sampleId]}
	publishDir "${barcodesOutputDir}", mode: 'move', overwrite: true

	when: params.link && !params.dryrun

	input:
	set sampleId, user from sampleAndUserLinkFastq
	file multiQcMarker from multiQcCompletionMarker

	output:
	file("${sampleId}_${user}_Link_Complete.txt") into linkCompletionMarkers

	script:
	"""
	mkdir -p /mnt/galaxy/tmp/users/${user}/${run_id}_shep && chmod 777 /mnt/galaxy/tmp/users/${user}/${run_id}_shep
	ln -sf ${barcodesOutputDir}/combined_fastq_tmp/${sampleId}*.fastq.gz /mnt/galaxy/tmp/users/${user}/${run_id}_shep
	touch ${sampleId}_${user}_Link_Complete.txt
	"""

}


process calculateBarcodeFrequency {
	tag { 'Calculating unassigned barcode frequency'}
	publishDir "${barcodesOutputDir}", mode: 'move', overwrite: true

	when: !params.dryrun

	input:
	file linkCompletionMarker from linkCompletionMarkers.toList()

	output:
	file('unassigned_barcode_frequency.txt') into barcode_distribution

	script:
	"""
	# Do not make the assumption there will always be tile 1101 and lane 1, particularly not the lane.
	# For instance, in NextSeq tile patterns are more like 11101 as opposed to 1101 for novaseq since there are a different count of lanes, surface, swath, tile in NextSeq
	zcat ${barcodesOutputDir}/fastq/*/L*_unassigned.barcode_*.fastq.gz |head -n 100000|sed -n '2~4p'|sort|uniq -c|sort -rn|head>unassigned_barcode_frequency.txt
	"""
}


process prepareEmailMessage {
	tag { 'Preparing email messages'}
	publishDir "${barcodesOutputDir}", mode: 'move', overwrite: true

	input:
	file barcodeDistribution from barcode_distribution

	output:
	file('message.txt') into email_message

	script:
	"""
	printf "Run ${experimentName} ( ${flowcell}) processing is complete. Output files have been linked to your user folder.\n\n \
	Output files have been linked to your user folder /mnt/galaxy/tmp/users/${user}/${run_id}_shep.\n\n \
	Barcode Summary from Lane 1 and the 10 most frequent barcodes observed in unassigned fastq file (first 100k reads) are attached.\n\n \
  Please remember to clean the instrument.">message.txt
	"""

}


process notifyUsersNF {

	tag { user_email }
	publishDir "${barcodesOutputDir}", mode: 'move', overwrite: true

	when: !params.nomail && params.mail_function.toLowerCase()=="nextflow" && !params.dryrun

	input:
	val user_email from user_emails_nextflow
	file barcode_freq from barcode_distribution
	file message from email_message
	file kronaHtml from kronaReport

	script:
	sendMail {
					from 'seq-shepherd@neb.com'
					to user_email
					cc 'ngebremedhin@neb.com'
          subject "Run Complete: ${experimentName} (${run_id})"
          attach "${barcodesOutputDir}/unassigned_barcode_frequency.txt", fileName: 'unassigned_barcode_frequency.txt'
          attach "${barcodesOutputDir}/barcode_summary.txt", fileName: 'barcode_summary.txt'
					attach "${barcodesOutputDir}/fastqc/multiqc_report.html", fileName: 'multiqc_report.html'
					attach "${barcodesOutputDir}/kraken/krona_report.html", fileName: 'krona_report.html'
          body """
					Run ${experimentName} ( ${flowcell}) processing is complete. Output files have been linked to your user folder.

					Output files have been linked to your user folder /mnt/galaxy/tmp/users/${user}/${run_id}_shep.

					Barcode Summary from Lane 1 and the 10 most frequent barcodes observed in unassigned fastq file (first 100k reads) are attached.

				  Please remember to clean the instrument.
          """
  }
	"""
	touch 'pipeline_complete.txt' # Marks the completion of the pipeline process
	"""
}

process notifyUsersRuby {
	tag { user_email }

	when: !params.nomail && params.mail_function.toLowerCase()=="ruby" && !params.dryrun

	input:
	val user_email from user_emails_ruby
	file barcode_freq from barcode_distribution
	file message from email_message
	file kronaHtml from kronaReport

	shell:
	"""
	ruby /mnt/home/ngebremedhin/novaseq_pipe/hts-scripts/nextflow-novaseq/mail.rb -r ${user_email} \
	-s  "Run Complete: ${experimentName} (${run_id})" -b ${barcodesOutputDir}/message.txt -f "${barcodesOutputDir}/unassigned_barcode_frequency.txt,\
	 ${barcodesOutputDir}/barcode_summary.txt, ${barcodesOutputDir}/fastqc/multiqc_report.html, ${barcodesOutputDir}/kraken/krona_report.html"
  touch 'pipeline_complete.txt' # Marks the completion of the pipeline processs
	"""
}

// End of Processes #############################################################################################################################


//++++++++++++++++++++++++++++++++++++++++++++++ Functions +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def cleanSampleId(sampleId) {
	return (sampleId =~ /^L._/).replaceFirst("")
}

def getReadsToTransfer(num_reads) {
    def readsToTransfer
    switch (num_reads) {
        case 1:
            readsToTransfer = [1]
            break
        case 2:
            readsToTransfer = [1]
            break
        case 3:
            readsToTransfer = [1,2] //The second read could be the third to be sequenced but still it is the second of two/paired end reads
            break
        case 4:
            readsToTransfer = [1,2]
            break
    }
    return readsToTransfer
}

def getRunDirectory(filePath) {
	x = filePath.tokenize("/")
	x.pop()
	return '/' +  x.join("/")

}

def createOutputDir(filePath) {

	def barcodesDir = new File(getRunDirectory(filePath) + "/Nextflow/Barcodes")
	if(!barcodesDir.exists()) {
		barcodesDir.mkdirs()
	} else {
		// Nextflow must have run before, so rename the outer directory
		// to-do: clean-up repetitive use of getRunDirectory
			new File(getRunDirectory(filePath) + "/Nextflow").renameTo(getRunDirectory(filePath) + "/Nextflow" + '_' + randomGenerator(7)) // rename the original outer 'Nextflow' directory
			barcodesDir.mkdirs()
	}
	return barcodesDir

}

def getTileXmlFile(filePath) {
	// A function that is meant to identify the xml file that contains list of tiles. Tiles are listed in RunInfo.xml in the
	// case of Novaseq and NextSeq while in config.xml in the case of MiSeqs.
	x = filePath.tokenize("/")
	x.pop()
	run_directory = x.join("/")


	runInfoXml = "/" + run_directory + "/RunInfo.xml"
	configXml = "/" + run_directory + "/Data/Intensities/BaseCalls/config.xml" // For Miseqs the tile info is in config.xml file buried under the Basecalls dir.

	if ( file(configXml).exists() ) { return [configXml, 'miseq'] } // Although both Novaseq and miseq have RunInfo.xml, this would suffice to tell which
	if ( file(runInfoXml).exists() ) { return [runInfoXml, 'novaseq'] }

	exit 1, "Missing the expected run xml file: filePath"
}

def parseRunXml(run_directory) {
	xmlFilePath = "/" + run_directory + "/RunInfo.xml" // The content varies but both Novaseq and Miseq have RunInfo.xml.
	runInfo = new XmlParser().parse(new File(xmlFilePath))
	num_reads = runInfo.Run.Reads.Read.size()
	run_id = runInfo.Run.@'Id'[0]
	flowcell = runInfo.Run.Flowcell.text()
	machine_name = runInfo.Run.Instrument.text()
	lanecount = runInfo.Run.FlowcellLayout.@'LaneCount'[0]

	read_structure = ""
	runInfo.Run.Reads.Read.each { r ->
		if ('Y'.equalsIgnoreCase(r.@'IsIndexedRead')) {
			read_structure += r.@'NumCycles' + "B"
		} else {
			read_structure += r.@'NumCycles' + "T"
		}

	}
	return [num_reads, run_id, flowcell, machine_name, lanecount, read_structure]
}

def retrieveListOfTiles(xmlFilePath, instrument, run_directory) {
	runInfo = new XmlParser().parse(new File(xmlFilePath))
	if ('novaseq'.equalsIgnoreCase(instrument)) {
		tiles = runInfo.Run.FlowcellLayout.TileSet.Tiles.Tile
	} else if('miseq'.equalsIgnoreCase(instrument)) {
		tiles = runInfo.Run.TileSelection.Lane.Tile
	} else {
		// Should think about how to handle other instrument types. It could be a NextSeq or Firefly or a new instrument
	}

	def preProcessDir = new File(run_directory + "/Nextflow/preProcess")
	def tilesFile = new File("${preProcessDir}/tilesList.csv")

	if(!preProcessDir.exists()) { preProcessDir.mkdirs()}

	tiles.each {
		r ->
		tilesFile << r.text() + "\n"
	}

	return [tilesFile, preProcessDir]
}

def parseSampleSheet(sampleSheetPath) {
  // Parse samplesheet.csv and return the experiment name and number of header row to skip
  header = 0
  experiment = 'NA'
  sampleSheetFile = new File(sampleSheetPath)
  sampleSheetFile.eachLine { line, lineNumber ->
    fields = line.split(',')
    //Randomly check rows for column names to demarcate end of header.
    if (fields.length>5 && fields[0].trim().equalsIgnoreCase('Sample_ID') && fields[1].trim().equalsIgnoreCase('Sample_Name') && fields[2].trim().equalsIgnoreCase('Sample_Plate') && fields[4].trim().equalsIgnoreCase('I7_Index_ID') ) {
        header = (lineNumber.toInteger() - 1)
    }
    if(fields.length>=2 && line.split(',')[0].trim().equalsIgnoreCase('Experiment Name')) {
      experiment = line.split(',')[1]
    }
  }
  return [header, experiment]
}

def randomGenerator(n) {
	alphanumeric_set = (('A'..'Z')+('0'..'9')).join()
	return new Random().with {
    (1..n).collect { alphanumeric_set[ nextInt( alphanumeric_set.length() ) ] }.join()
  }
}
//End of Functions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
