//
// ExternalSort.groovy
//
// Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
//

#!/usr/bin/env groovy

import groovy.transform.Sortable
import groovy.transform.ToString

import java.util.zip.Deflater
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

/**
 * ExternalSort.groovy
 * 
 * Sorts reads of a fastq or fasta file according to their ID.
 * If the sorting cannot be done entirely in memory
 * a divide and conquer approach is followed.
 * Input data is split in n parts(depending the size of data and available memory),  
 * Each part is been sorted in memory.
 * Each sorted part is saved in a temporary file.
 * The sorted temporary files are merged to a final sorted output file
 * 
 * Input fastq reads should be of 4 lines each (header,sequence,+,qualities)
 * Input fasta reads should be of 2 lines each (header,sequence)
 * 
 * @author Anestis Gkanogiannis
 * @email ganoyan@gmail.com
 * @version 1.0
 * @date 26/05/2016
 *
 **/
class ExternalSort {
	static final def DEFAULTMAXTMPFILES = 1024
	static final def DEFAULTTMPSTORE = '/tmp'
	
	static def defaultcomparator = new Comparator<Read>() {
		int compare(Read r1, Read r2) {
			return r1.compareTo(r2)
		}
	}
	
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
		def externalsort = new ExternalSort()
		externalsort.execute(args)
	}
	
	void execute(def args) {
		def cli = new CliBuilder(usage:'groovy ExternalSort.groovy [options] [input] [output]', header:'Options:')
		cli.with {
			h longOpt: 'help', 'Show usage information'
			f longOpt: 'fasta', 'Input is in fasta format otherwise is fastq, default=false'
			z longOpt: 'gzip', 'Use compression for the temporary files, default=false'
			c longOpt: 'clean', 'Remove comments from header (#abcdefg), default=false'
			t longOpt: 'maxtmpfiles', args:1, 'Maximum number of intermediate temporary files, default=1024'
			s longOpt: 'tmpstore', args:1, 'Where to save the intermediate temporary files, default=/tmp'
		}
	
		def fasta = false
		def usegzip = false
		def clean = false
		def maxtmpfiles = DEFAULTMAXTMPFILES
		def tempfilestore = new File(DEFAULTTMPSTORE)
		def inputfile, outputfile
		
		def options = cli.parse(args)
		
		if(options.h) {
			cli.usage()
			return
		}
		
		if (options.f) fasta = true
		if (options.z) usegzip = true
		if (options.c) clean = true
		if (options.t) maxtmpfiles = Integer.parseInt(options.t)
		if (options.s) tempfilestore = new File(options.s)
		def extraArguments = options.arguments()
		if (extraArguments) {
			inputfile = new File(extraArguments[0])
			if (extraArguments.size() > 1) 
				outputfile = new File(extraArguments[1])
		}
		
		if (inputfile==null || outputfile==null) {
			println  "Please provide input and output file names."
			cli.usage()
			return
		}
		
		def isgzipped = isGZipped(inputfile)
		
		println 'Available memory (bytes) :' + estimateAvailableMemory()
		println 'Input is ' + (isgzipped?'':'not ') + 'gzipped.'
		
		def br
		if(isgzipped)
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputfile), 2048)))
		else
			br = new BufferedReader(new FileReader(inputfile))
		def rbr = new ReadBufferedReader(br: br, clean: clean, fasta: fasta)
		def datalength
		if(isgzipped)
			datalength = 5 * inputfile.length()
		else
			datalength = inputfile.length()
		def sortedtmpfiles = sort(rbr, datalength, defaultcomparator, estimateAvailableMemory(), usegzip, maxtmpfiles, tempfilestore)
		println 'Created ' + sortedtmpfiles.size() + ' tmp files.'
		def readcounter = merge(sortedtmpfiles, outputfile, defaultcomparator, fasta, usegzip)
		println 'Sorted ' + readcounter + ' reads.'
	}
	
	def estimateAvailableMemory() {
		System.gc()
		return Runtime.getRuntime().freeMemory()
	}

	def estimateBestSizeOfBlocks(long sizeoffile, int maxtmpfiles, long maxmemory) {
		def blocksize = sizeoffile / maxtmpfiles + (sizeoffile % maxtmpfiles == 0 ? 0 : 1)
		if (blocksize < maxmemory/2) blocksize = maxmemory/2
		return (long)blocksize
	}
	
	def sort(ReadBufferedReader rbr, long datalength, Comparator cmp, long maxMemory, 
		boolean usegzip, int maxtmpfiles, File tempFileStore) {
		def sortedtmpfiles = []
		def blocksize = estimateBestSizeOfBlocks(datalength, maxtmpfiles, maxMemory)
		
		try {
			def tmplist = []
			def read = new Read("", "")
			try {
				while (read != null) {
					def currentblocksize = 0L
					while ((currentblocksize < blocksize) && ((read = rbr.getRead()) != null)) {
						tmplist.add(read)
						currentblocksize += (long)StringSizeEstimator.estimatedSizeOf(read)
					}
					sortedtmpfiles.add(sortAndSave(tmplist, cmp, tempFileStore, usegzip, sortedtmpfiles.size()+1))
					tmplist.clear()
				}
			} 
			catch (all) {
				if (tmplist.size() > 0) {
					sortedtmpfiles.add(sortAndSave(tmplist, cmp, tempFileStore, usegzip, sortedtmpfiles.size()+1))
					tmplist.clear()
				}
			}
		} 
		finally {
			rbr.br.close()
		}
		return sortedtmpfiles
	}
		
	def sortAndSave(List tmplist, Comparator cmp, File tempfilestore, boolean usegzip, int tmpfilescounter) {
		//Java 8
		//tmplist.sort(cmp)
		//Java 7
		Collections.sort(tmplist, cmp)
		def newtmpfile = File.createTempFile('ExternalSort', 'tmpFile', tempfilestore)
		println 'Saving tmpfile:'+ tmpfilescounter + ' (' + newtmpfile.getCanonicalPath() + ')'
		newtmpfile.deleteOnExit()
		def out = new FileOutputStream(newtmpfile)
		def ZIPBUFFERSIZE = 2048
		if (usegzip)
			out = new GZIPOutputStream(out, ZIPBUFFERSIZE) {{this.def.setLevel(Deflater.BEST_SPEED)}}
		def bw = new BufferedWriter(new OutputStreamWriter(out))
		try {
			for (Read r : tmplist) 
				bw.write(r.header+'\n'+r.data+'\n')
		} 
		finally {
			bw.close()
		}
		return newtmpfile
	}
	
	def merge(List sortedfiles, File outputfile, Comparator cmp, boolean fasta, boolean usegzip) {
		def readbuffers = []
		for (File f : sortedfiles) {
			def ZIPBUFFERSIZE = 2048
			def is = new FileInputStream(f)
			def br
			if(usegzip)
				br = new BufferedReader(new InputStreamReader(new GZIPInputStream(is, ZIPBUFFERSIZE)))
			else
				br = new BufferedReader(new InputStreamReader(is))
			readbuffers.add(new ReadBuffer(br, fasta))
		}
		def bw = new BufferedWriter(new FileWriter(outputfile))
		def readcounter = mergeSortedFiles(bw, cmp, readbuffers)
		for (File f : sortedfiles) f.delete()
		return readcounter
	}
	
	def mergeSortedFiles(BufferedWriter bw, Comparator cmp, List readbuffers) {
		PriorityQueue<ReadBuffer> pq = new PriorityQueue<ReadBuffer>(11, 
			new Comparator<ReadBuffer>() {
				int compare(ReadBuffer i, ReadBuffer j) {
					return cmp.compare(i.peek(), j.peek())
				}
			})
		for (ReadBuffer readbuffer : readbuffers)
			if (!readbuffer.empty())
				pq.add(readbuffer)
		def readcounter = 0
		try {
			while (pq.size() > 0) {
				ReadBuffer readbuffer = pq.poll()
				Read r = readbuffer.pop()
				bw.write(r.header+'\n'+r.data+'\n')
				++readcounter;
				if (readbuffer.empty())
					readbuffer.br.close()
				else
					pq.add(readbuffer)
			}
		}
		finally {
			bw.close()
			for (ReadBuffer rb : pq)
				rb.br.close()
		}
		return readcounter
	}
}

@Sortable(includes=['id'])
@ToString
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
}

class ReadBufferedReader {
	BufferedReader br
	def clean = false
	def fasta = false

	def getRead() {
		def read
		String line
		while((line=br.readLine())!=null) {
			if(!fasta && line.startsWith('@')) {
				if(clean && line.contains('#'))
					line = line.replaceAll('#.*?(?=\\/|$)','')
				read = new Read(line, br.readLine()+'\n'+br.readLine()+'\n'+br.readLine())
				break
			}
			else if(fasta && line.startsWith('>')) {
				if(clean && line.contains('#'))
					line = line.replaceAll('#.*?(?=\\/|$)','')
				read = new Read(line, br.readLine())
				break
			}
		}
		return read
	}
}

class ReadBuffer {
	BufferedReader br
	Read cache
	def fasta = false

	public ReadBuffer(BufferedReader br, boolean fasta) {
		this.br = br
		this.fasta = fasta
		reload()
	}
	
	def empty() {
		return cache == null;
	}
	
	def peek() {
		return cache;
	}
	
	def pop() {
		def answer = new Read(peek())
		reload()
		return answer
	}
	
	void reload() {
		cache = null
		def line
		while((line=br.readLine())!=null) {
			if(!fasta && line.startsWith('@')) {
				cache = new Read(line, br.readLine()+'\n'+br.readLine()+'\n'+br.readLine())
				break
			}
			else if(fasta && line.startsWith('>')) {
				cache = new Read(line, br.readLine())
				break
			}
		}
	}

}

class StringSizeEstimator {
	static int OBJ_OVERHEAD
	
	static {
		def IS_64_BIT_JVM = true
		def arch = System.getProperty("sun.arch.data.model")
		if (arch!=null && arch.indexOf("32")!=-1) IS_64_BIT_JVM = false
		def OBJ_HEADER = IS_64_BIT_JVM ? 16 : 8
		def ARR_HEADER = IS_64_BIT_JVM ? 24 : 12
		def OBJ_REF = IS_64_BIT_JVM ? 8 : 4
		def INT_FIELDS = 12
		OBJ_OVERHEAD = OBJ_HEADER + INT_FIELDS + OBJ_REF + ARR_HEADER
	}
	
	static def estimatedSizeOf(String s) {
		return (s.length() * 2) + OBJ_OVERHEAD
	}
		
	static def estimatedSizeOf(Read r) {
		return ((r.id.length()+r.header.length()+r.data.length()) * 2) + OBJ_OVERHEAD
	}
}
