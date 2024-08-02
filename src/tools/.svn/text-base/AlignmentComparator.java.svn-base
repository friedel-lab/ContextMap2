package tools;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.util.Date;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NavigableMap;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

public class AlignmentComparator<E extends Comparable<E>> {

	public static void main(String args[]) {
		try {
			if(args[0].equals("getoverlap")) {
				String cm_readids_path = null;
				String bwt_readids_path = null;
				String outputfile_path = null;
				for(int i = 1; i < args.length; i++) {
					if(args[i].equals("-cmids")) {
						cm_readids_path = args[++i];
						continue;
					}
					if(args[i].equals("-bwtids")) {
						bwt_readids_path = args[++i];
						continue;
					}
					if(args[i].equals("-o")) {
						outputfile_path = args[++i];
						continue;
					}
				}
				getOverlap(cm_readids_path,bwt_readids_path,outputfile_path);
			}
			else if(args[0].equals("getsequences")) {
				String idsPath = null;
				String fastaPath = null;
				String outputfile_path = null;
				for(int i = 1; i < args.length; i++) {
					if(args[i].equals("-ids")) {
						idsPath = args[++i];
						continue;
					}
					if(args[i].equals("-fasta")) {
						fastaPath = args[++i];
						continue;
					}
					if(args[i].equals("-o")) {
						outputfile_path = args[++i];
						continue;
					}
				}
				getReadSequences(fastaPath, idsPath, outputfile_path);
			}
			else if(args[0].equals("filtersam")) {
				String samPath = null;
				String fastaPath = null;
				String idsPath = null;
				String outputPath = null;
				boolean positiveFilter = true;
				for(int i = 1; i < args.length; i++) {
					if(args[i].equals("-sam")) {
						samPath = args[++i];
						continue;
					}
					if(args[i].equals("-fasta")) {
						fastaPath = args[++i];
						continue;
					}
					if(args[i].equals("-ids")) {
						idsPath = args[++i];
						continue;
					}
					if(args[i].equals("-o")) {
						outputPath = args[++i];
						continue;
					}
					if(args[i].equals("--negative")) {
						positiveFilter = false;
						continue;
					}
				}
				filterSam(samPath,fastaPath,idsPath, outputPath,positiveFilter);
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void filterSam(String samPath, String fastaPath, String idsPath, String outputPath, boolean positiveFilter) {
		try {
			BufferedReader br; 
			HashSet<String> ids = new HashSet<String>();
			String currentLine;
			
			if(fastaPath != null) {
				br	= new BufferedReader(new FileReader(new File(fastaPath)));
				while((currentLine = br.readLine()) != null) {
					if(currentLine.charAt(0) == '>')
						ids.add(currentLine.substring(1));
				}
				br.close();
			}
			else {
				br	= new BufferedReader(new FileReader(new File(idsPath)));
				while((currentLine = br.readLine()) != null) {
					ids.add(currentLine);
				}
				br.close();
			}
			
			br = new BufferedReader(new FileReader(new File(samPath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				if(positiveFilter && ids.contains(splittedLine[0]))
					pw.println(currentLine);
				else if(!positiveFilter && !ids.contains(splittedLine[0]))
					pw.println(currentLine);
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void getOverlap(String cm_readids_path, String bwt_readids_path, String outputfile_path) {
		try {
			HashSet<String> cm_readids = new HashSet<String>();
			BufferedReader br = new BufferedReader(new FileReader(new File(cm_readids_path)));
			String currentLine;
			
			while((currentLine = br.readLine()) != null) {
				cm_readids.add(currentLine);
				
			}
		    br.close();
		    int cm_reads = cm_readids.size();
		    int overlappingReadIds = 0;
		    int bwtReadIds = 0;
		    
		    br = new BufferedReader(new FileReader(new File(bwt_readids_path)));
		    while((currentLine = br.readLine()) != null) {
				if(cm_readids.contains(currentLine)) {
					overlappingReadIds++;
					cm_readids.remove(currentLine);
				}
				
				bwtReadIds++;
			}
		    br.close();
			
		    
		    PrintWriter pw = new PrintWriter(new FileWriter(new File(outputfile_path)));
		    for(String readid : cm_readids) {
		    	pw.println(readid);
		    }
		    pw.close();
		    
		    System.out.println("cm read ids:\t" + cm_reads);
		    System.out.println("bwt read ids:\t" + bwtReadIds);
		    System.out.println("overlap:\t" + overlappingReadIds);
		    
		    
		}
		catch(Exception e ) {
			e.printStackTrace();
		}
	}
	
	private static void getReadSequences(String fastaFilePath, String readIdsPath, String outputFilePath) {
		try {
			BufferedReader idsReader = new BufferedReader(new FileReader(new File(readIdsPath)));
			String currentLine;
			HashSet<String> ids = new HashSet<String>();
			while((currentLine = idsReader.readLine()) != null) {
				ids.add(currentLine);
			}
			idsReader.close();
			
			BufferedReader fastaReader = new BufferedReader(new FileReader(new File(fastaFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			while((currentLine = fastaReader.readLine()) != null) {
				if(currentLine.charAt(0) == '>' && ids.contains(currentLine.substring(1))) {
					pw.println(currentLine);
					pw.println(fastaReader.readLine());
				}
			}
			fastaReader.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
}
