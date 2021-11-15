#!/usr/bin/env groovy

import java.io.BufferedReader;
import java.util.zip.Deflater
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

/**
 * SynchronizePairReads.groovy
 * 
 * Synchronizes paired reads.
 * Pair reads have an ID of the form /1 and /2
 * Reads without a pair are written to _single.synced.fasta(q)
 * Forward reads are written to _1.synced.fasta(q)
 * Reverse reads are written to _2.synced.fasta(q)
 * 
 * Input reads should be sorted according to their ID.
 * Input can be fasta , fastq, gzipped or not.
 *
 * Input fastq reads should be of 4 lines each (header,sequence,+,qualities)
 * Input fasta reads should be of 2 lines each (header,sequence)
 *
 * @author Anestis Gkanogiannis
 * @email ganoyan@gmail.com
 * @version 1.0
 * @date 27/05/2016
 *
 **/
class SynchronizePairReads {
	
	static def isGZipped(File f) {
		def br = new BufferedInputStream(new FileInputStream(f))
		br.mark(2)
		def magic = 0
		try {
			magic = br.read() & 0xff | ((br.read() << 8) & 0xff00)
			br.reset()
			br.close()
		}
		catch (all) {
			return false
		}
		return magic == GZIPInputStream.GZIP_MAGIC
	}
	
	static main(def args) {
		def synchronizepaireads = new SynchronizePairReads()
		synchronizepaireads.execute(args)
	}
	
	void execute(def args) {
		def cli = new CliBuilder(usage:'groovy SynchronizePairReads.groovy [options] -f input_1 -r input_2', header:'Options:')
		cli.with {
			h longOpt: 'help', 'Show usage information'
			q longOpt: 'fastq', 'Input is in fastq format otherwise is fasta, default=false'
			z longOpt: 'gzip', 'Use compression for the output files, default=false'
			c longOpt: 'clean', 'Remove comments from header (#abcdefg), default=false'
			o longOpt: 'prefix', args:1, 'Output files folder'
			p longOpt: 'prefix', args:1, 'Prefix of output files'
			f longOpt: 'forward', args:1, 'Forward reads input file'
			r longOpt: 'reverse', args:1, 'Reverse reads input file'
		}
		
		def fastq = false
		def usegzip = false
		def clean = false
		def output
		def prefix
		def inputfile1, inputfile2
		
		def options = cli.parse(args)
		
		if(options.h) {
			cli.usage()
			return
		}
		
		if (options.q) fastq = true
		if (options.z) usegzip = true
		if (options.c) clean = true
		if (options.o) output = options.o
		if (options.p) prefix = options.p
		if (options.f) inputfile1 = new File(options.f)
		if (options.r) inputfile2 = new File(options.r)
		
		if (inputfile1==null || inputfile2==null) {
			println  "Please provide input_1 and input_2 file names."
			cli.usage()
			return
		}
		
		println 'input_1 is ' + (isGZipped(inputfile1)?'':'not ') + 'gzipped.'
		println 'input_2 is ' + (isGZipped(inputfile2)?'':'not ') + 'gzipped.'
		
		def br1, br2
		if(isGZipped(inputfile1))
			br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputfile1), 2048)))
		else
			br1 = new BufferedReader(new FileReader(inputfile1))
		def rbr1 = new ReadBufferedReader(br: br1, clean: clean, fastq: fastq)
		if(isGZipped(inputfile2))
			br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputfile2), 2048)))
		else
			br2 = new BufferedReader(new FileReader(inputfile2))
		def rbr2 = new ReadBufferedReader(br: br2, clean: clean, fastq: fastq)
		
		def suffix = fastq?'.synced.fastq':'.synced.fasta'
		def bw1 = new BufferedWriter(new FileWriter(new File(output, prefix+'_1'+suffix)))
		def bw2 = new BufferedWriter(new FileWriter(new File(output, prefix+'_2'+suffix)))
		def bwsingle = new BufferedWriter(new FileWriter(new File(output, prefix+'_single'+suffix)))
		
		try {
			def read1 = rbr1.getRead()
			def read2 = rbr2.getRead()
			while (read1 != null && read2 != null) {
				if(read1.compareTo(read2)>0) {
					bwsingle.write(read2.header+'\n'+read2.data+'\n')
					read2 = rbr2.getRead()
				}
				else if(read1.compareTo(read2)<0) {
					bwsingle.write(read1.header+'\n'+read1.data+'\n')
					read1 = rbr1.getRead()
				}
				else {
					bw1.write(read1.header+'\n'+read1.data+'\n')
					bw2.write(read2.header+'\n'+read2.data+'\n')
					read1 = rbr1.getRead()
					read2 = rbr2.getRead()
				}
			}
			if(read1!=null){
				while(true) {
					bwsingle.write(read1.header+'\n'+read1.data+'\n')
					if((read1=rbr1.getRead())==null) 
						break
				}
			}
			if(read2!=null){
				while(true) {
					bwsingle.write(read2.header+'\n'+read2.data+'\n')
					if((read2=rbr2.getRead())==null)
						break
				}
			}
		}
		finally {
			rbr1.br.close()
			rbr2.br.close()
			bw1.close()
			bw2.close()
			bwsingle.close()
		}
	}
}

class Read {
	String id
	String header
	String data
	
	Read(String header, String data) {
		this.header = header.toUpperCase()
		this.data = data
		if(this.header != null) {
			try {
				this.id = this.header.split("\\s+")[0].substring(1).toUpperCase().trim()
			}
			catch(e) {
				this.id = this.header.toUpperCase()
			}
		}
	}

	Read(Read other) {
		if(other!=null) {
			this.id = other.id.toString()
			this.header = other.header.toString().toUpperCase()
			this.data = other.data.toString()
		}
	}
	
	def compareTo(Read other) {
		def idthis = this.id.replaceAll('\\/.*$','')
		def idother = other.id.replaceAll('\\/.*$','')
		return idthis.compareTo(idother)
	}
}

class ReadBufferedReader {
	BufferedReader br
	def clean = false
	def fastq = false

	def getRead() {
		def read
		String line
		while((line=br.readLine())!=null) {
			if(fastq && line.startsWith('@')) {
				if(clean && line.contains('#'))
					line = line.replaceAll('#.*?(?=\\/|$)','')
				read = new Read(line, br.readLine()+'\n'+br.readLine()+'\n'+br.readLine())
				break
			}
			else if(!fastq && line.startsWith('>')) {
				if(clean && line.contains('#'))
					line = line.replaceAll('#.*?(?=\\/|$)','')
				read = new Read(line, br.readLine())
				break
			}
		}
		return read
	}
}