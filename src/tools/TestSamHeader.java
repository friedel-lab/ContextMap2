package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import alignment.SamStreamConsumer;
import alignment.StreamConsumerSimple;
import alignment.StreamType;

/**
 * Functions to check, if the reference genome is available for all names stored in the index file.
 *
 */
public class TestSamHeader {
	
	private static final String[] FASTA_ENDINGS = {"fa", "fasta"};
	private static final String SAM_INDEX_REGEX = "^@SQ\tSN:([^\t]+)\tLN:([0-9]+).*";
	private static final String SAM_PROGRAM_TAG = "@PG";
	private static final String NAME_OF_DUMMY_FASTA_FILE = "dummyFastaFile.fa";
	private static final String NAME_OF_SAM_TEST_INDEX_FILE = "testIndex_dummyRead.sam";
	private static final String DUMMY_READ = ">dummyRead1" + System.lineSeparator() + "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA";
	private static final String FASTA_HEADER_START = ">";

	/**
	 * Test method
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		ArrayList<String> reference = getNameOfGenomeFiles("/mnt/biostor1/Data/Databases/GENOMES/Homo/simplefa");
		//ArrayList<String> index = getReferenceSequenceNameFromHeader("/home/proj/biosoft/software/bowtie-1.0.0/bowtie", "/mnt/biostor1/kluge/bowtie_index_test/bowtie1", "/tmp/headerTest", null);
		//ArrayList<String> index = getReferenceSequenceNameFromHeader("/home/proj/biosoft/software/bowtie2-2.1.0/bowtie2", "/mnt/biostor1/kluge/bowtie_index_test/bowtie2", "/tmp/headerTest", null);
		ArrayList<String> index = getReferenceSequenceNameFromHeader("/home/proj/biosoft/software/bwa-0.7.9a/bwa", "/mnt/biostor1/kluge/test.fa", "/tmp/headerTest", null);
		ArrayList<String> diff = compareIndexAndReference(index, reference);
		
		System.out.println("Names in reference:");
		for(String r : reference) {
			System.out.println(r);
		}
		System.out.println("______________________");
		System.out.println("Names in index:");
		for(String i : index) {
			System.out.println(i);
		}
		System.out.println("______________________");
		if(diff.size() == 0) {
			System.out.println("All reference sequences were found (" + diff.size() + ").");
		}
		else {
			System.out.println(diff.size()+ " reference sequences are missing:");
		}
		for(String d : diff) {
			System.out.println(d);
		}
	}
	
	/** bowtie1 **/
	/*private static String getSamHeaderCommand(String alignerBinPath, String genomeIndexBasePath, String dummyReadPath) {
		return String.format("%s -f --sam %s %s", alignerBinPath, genomeIndexBasePath, dummyReadPath);
	}*/
	
	/** bowtie2 **/
	/*private static String getSamHeaderCommand(String alignerBinPath, String genomeIndexBasePath, String dummyReadPath) {
		return String.format("%s -f %s %s", alignerBinPath, genomeIndexBasePath, dummyReadPath);
	}*/
	
	/** bwa **/
	private static String getSamHeaderCommand(String alignerBinPath, String genomeIndexBasePath, String dummyReadPath) {
		return String.format("%s mem %s %s", alignerBinPath, genomeIndexBasePath, dummyReadPath);
	}
	
	public ArrayList<String> checkIndex(String alignerBinPath, String genomeIndexBasePath, String referenceDir, String tmpFolder, ReadAligner aligner) throws Exception {
		ArrayList<String> reference = getNameOfGenomeFiles(referenceDir);
		ArrayList<String> index = getReferenceSequenceNameFromHeader(alignerBinPath, genomeIndexBasePath, tmpFolder, aligner);
		ArrayList<String> diff = compareIndexAndReference(index, reference);
		
		return diff;
	}

	/**
	 * Try to get a SAM header from aligner tool by aligning of a dummy read file
	 * @param alignerBinPath
	 * @param genomeIndexBasePath
	 * @param tmpFolder
	 * @param aligner
	 * @return
	 * @throws IOException
	 */
	public static File getSAMHeader(String alignerBinPath, String genomeIndexBasePath, String tmpFolder, ReadAligner aligner) throws IOException {
		// create tmp folder, if not there
		File tmp = new File(tmpFolder);
		// create dummy read in tmp folder
		if((tmp.exists() || tmp.mkdirs())) {	
			File dummyRead = createDummyReadFile(tmpFolder);
			if(dummyRead != null) {
				File samTestIndexOutput = new File(tmpFolder + File.separator + NAME_OF_SAM_TEST_INDEX_FILE);
				samTestIndexOutput.delete();
				// try to get a dummy SAM index
				try {
					Process process;
					int returnValue;
					SamStreamConsumer errorOut, stdOut;	
						
					// call tool
					process = Runtime.getRuntime().exec(getSamHeaderCommand(alignerBinPath, genomeIndexBasePath, dummyRead.getAbsolutePath()));
					// TODO: call this one!
					//process = Runtime.getRuntime().exec(aligner.getSamHeaderCommand(alignerBinPath, genomeIndexBasePath, dummyRead.getAbsolutePath()));
					errorOut = new StreamConsumerSimple(process.getErrorStream(), StreamType.ERROR, null);
					stdOut = new StreamConsumerSimple(process.getInputStream(), StreamType.STDOUT, samTestIndexOutput.getAbsolutePath());
				
					errorOut.start();
					stdOut.start();
					returnValue = process.waitFor();
				
					//notify stream threads that the process terminated
					errorOut.processTerminated();
					stdOut.processTerminated();
				
					synchronized(errorOut) {
						while(errorOut.isAlive()) {
							errorOut.wait();
						}
					}
					synchronized(stdOut) {
						while(stdOut.isAlive()) {
							stdOut.wait();
						}
					}
							
					process.getErrorStream().close();
					process.getOutputStream().close();
					process.getInputStream().close();
					process.destroy();
				
					if(returnValue != 0) {
						throw new Exception(String.format("Something went wrong while getting sam header, aligner tool return value: %s.",returnValue));
					}	
					else {
						// all should be ok
						if(samTestIndexOutput.exists() && samTestIndexOutput.canRead())
							return samTestIndexOutput;
					}
				}
				catch(Exception e) {
					e.printStackTrace();
				}
			}
		}
		else {
			System.out.println("Could not create tmp folder '"+ tmp.getAbsolutePath() +"'");
		}
		return null;
	}
	
	
	/**
	 * Create a dummy read file
	 * @param tmpFolder
	 * @return
	 * @throws IOException
	 */
	public static File createDummyReadFile(String tmpFolder) throws IOException {
		// create tmp folder, if not there
		File tmp = new File(tmpFolder);
		if(tmp.exists() || tmp.mkdirs()) {	
			// create dummy read in tmp folder
			File dummyRead = new File(tmpFolder + File.separator + NAME_OF_DUMMY_FASTA_FILE);
			dummyRead.delete();
			if(dummyRead.createNewFile()) {
				FileWriter fw = new FileWriter(dummyRead);
				fw.write(DUMMY_READ);
				fw.flush();
				fw.close();
			}
			return dummyRead;
		}
		return null;
	}
	

	/**
	 * Finds the index names from the SAM header
	 * @param alignerBinPath
	 * @param genomeIndexBasePath
	 * @param tmpFolder
	 * @param aligner
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<String> getReferenceSequenceNameFromHeader(String alignerBinPath, String genomeIndexBasePath, String tmpFolder, ReadAligner aligner) throws Exception {
		// call aligner in order to get the sam header
		File samHeader = getSAMHeader(alignerBinPath, genomeIndexBasePath, tmpFolder, aligner);
	
		if(samHeader != null) {
			ArrayList<String> files = new ArrayList<String>();
			Pattern p = Pattern.compile(SAM_INDEX_REGEX);
			// open the file
			BufferedReader br = new BufferedReader(new FileReader(samHeader));
			
			// find lines which are matching the regex pattern
			String line, name;
			@SuppressWarnings("unused")
			int length; // a more strict implementation might also want to test, if the length of both sequenes is equal
			
			while((line = br.readLine()) != null) {
				 Matcher m = p.matcher(line);
				 // pattern matched!
				 if(m.find()) {
					 name = m.group(1);
					 length = Integer.parseInt(m.group(2));
					 files.add(name);
				 }
				 // end search if @PG tag was found
				 else if(line.startsWith(SAM_PROGRAM_TAG)) {
					break;
				 }
			}
			// close the buffer and return the found names
			br.close();
			
			return files;
		}
		else {
			System.out.println("SamHeader was not created successfully...)");
		}
		return null;
	}
	
	
	/**
	 * Finds all fasta files in that folder and returns the header names of the fasta header 
	 * @param referenceDir
	 * @return
	 */
	public static ArrayList<String> getNameOfGenomeFiles(String referenceDir) {
		ArrayList<String> names = new ArrayList<String>();
		File refpath = new File(referenceDir);
		
		String tmp[];
		// check, if path is a folder
		if(refpath.isDirectory()) {
			// get all files in folder
			for(String name : refpath.list()) {
				// check, if file ending is valid
				tmp = name.split("\\.");
				for(String vEnding : FASTA_ENDINGS) {
					if(tmp[tmp.length-1].equals(vEnding)) {
						// read the fasta file and try to find headers
						try {
							File f = new File(referenceDir + File.separator + name);
							if(f.isFile() && f.canRead()) {
								BufferedReader br = new BufferedReader(new FileReader(f));
		
								String line, hname;
								while((line = br.readLine()) != null) {
									// test if line is fasta header
									if(line.startsWith(FASTA_HEADER_START)) {
										hname = line.replaceFirst(FASTA_HEADER_START + "\\s*", "").split("\\s")[0];
										names.add(hname);
										
										System.out.println("Genomedir: Filename: " + name + " --- headname: " + hname);
									}
								}
								// close the buffer and return the found names
								br.close();
							}
						}
						catch(IOException e) {
							e.printStackTrace();
						}
						// one valid ending is enough
						break;
					}
				}
			}
		}
		return names;
	}
	
	
	/**
	 * Compares two lists of names. It is check, if all entries in index are also stored in reference
	 * @param index
	 * @param reference
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static ArrayList<String> compareIndexAndReference(ArrayList<String> index, ArrayList<String> reference) throws IllegalArgumentException {
		// check, if input is valid
		if(index != null && reference != null) {
			ArrayList<String> missing = new ArrayList<String>();
			HashSet<String> ref = new HashSet<String>(reference);
			
			// check, if all what is in the index is also in the reference folder
			for(String i : index) {
				if(!ref.contains(i))
					missing.add(i);
			}
			return missing;
		}
		throw new IllegalArgumentException("None of the input parameters can be null.");
	}
}
