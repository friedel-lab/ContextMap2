package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;

import static java.nio.file.StandardCopyOption.*;

import java.nio.channels.FileChannel;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Date;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;



public class FileHandler {

	
	public FileHandler() {
		
	}
	
	public enum FileType {
		ZIP, GZ, TAR, TARGZ, FASTA, FASTQ, UNKNOWN
	}
	
	public enum ReadFormat {
		FASTA, FASTQ, UNKNOWN
	}
	
	public void checkInput(ArrayList<String> filePaths) {
		Date date;
		for(String readFilePath : filePaths) {
			if(!new File(readFilePath).isFile()) {
				date = new Date();
				System.err.println(String.format("[%s]\t%s is not a valid file path. Please check the help by calling ContextMap without any parameters",date.toLocaleString(), readFilePath));
				System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
				System.exit(1);
			}
			try {
				BufferedReader tmpReader = new BufferedReader(new FileReader(new File(readFilePath)));
				if(tmpReader.readLine() == null) {
					date = new Date();
					System.err.println(String.format("[%s]\tFound an empty read sequences file (%s). Please check the help by calling ContextMap without any parameters",date.toLocaleString(), readFilePath));
					System.err.println(String.format("[%s]\tAborting ContextMap run.",date.toLocaleString()));
					System.exit(1);
				}
				tmpReader.close();
			}
			catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
	}
	
	
	public void processInput(ArrayList<String> filePaths, String tmpOutDir, String outputPath) {
		try {
			
			if(!new File(tmpOutDir).exists())
				new File(tmpOutDir).mkdirs();
			
			FileType fileType;
			ReadFormat readFormat;
			File tmpOutputFile;
			File inputFile;
			ArrayList<String> tmpFilePaths = new ArrayList<String>();
			ArrayList<Integer> archiveFilesIndices = new ArrayList<Integer>();
			int numberOfFastaFiles = 0;
			//extracting archive files and converting fastq to fasta files
			for(String filePath : filePaths) {
				fileType = getFileType(filePath);
				tmpOutputFile = File.createTempFile("input", null, new File(tmpOutDir));
				inputFile = new File(filePath);
				tmpFilePaths.add(tmpOutputFile.getAbsolutePath());
				switch(fileType) {
					case ZIP:
						unzip(inputFile,tmpOutputFile);
						archiveFilesIndices.add(tmpFilePaths.size() - 1);
						break;
						
					case GZ:
						gunzip(inputFile,tmpOutputFile);
						archiveFilesIndices.add(tmpFilePaths.size() - 1);
						break;
						
					case TAR:
						//TODO
						throw new Exception("TAR archives are not supported yet.");
						//archiveFilesIndices.add(tmpFilePaths.size() - 1);
						//break;
						
					case TARGZ:
						//TODO 
						throw new Exception("TAR.GZ archives are not supported yet.");
						//archiveFilesIndices.add(tmpFilePaths.size() - 1);
						//break;
						
					case FASTQ:
						convertFastqToFasta(inputFile,tmpOutputFile);
						numberOfFastaFiles++;
						break;
						
					case FASTA:
						tmpFilePaths.remove(tmpFilePaths.size() - 1);
						tmpFilePaths.add(filePath);
						numberOfFastaFiles++;
						break;
						
					case UNKNOWN:
						throw new Exception(String.format("Unknown file extension for file: %s. Aborting ContextMap run.", filePath));
				}
			}
			
			
			//converting decompressed fastq files to fasta files
			
			for(int index : archiveFilesIndices) {
				readFormat = getReadFormat(tmpFilePaths.get(index));
				switch(readFormat) {
					case FASTA:
						numberOfFastaFiles++;
						break;
				
					case FASTQ:
						tmpOutputFile = File.createTempFile("input", null, new File(tmpOutDir));
						convertFastqToFasta(new File(tmpFilePaths.get(index)),tmpOutputFile);
						numberOfFastaFiles++;
						tmpFilePaths.set(index, tmpOutputFile.getAbsolutePath());
						break;
						
					case UNKNOWN:
						throw new Exception(String.format("Unknown file format in input archive found: %s. Aborting ContextMap run.", filePaths.get(index)));
				}
			}
			
			
			//all files are now in fasta format.
			//Currently only one or two input files are allowed: one input file -> single end mode . two input files -> paired end mode
			//code can be easily adapted to support more than two input files.
			if(tmpFilePaths.size() == numberOfFastaFiles) {
				//single end
				if(numberOfFastaFiles == 1) {
					//we generated a tmp fasta file and can simply move it to the final output path
					if(!tmpFilePaths.get(0).equals(filePaths.get(0))) {
						Files.move(FileSystems.getDefault().getPath(tmpFilePaths.get(0)),FileSystems.getDefault().getPath(outputPath),REPLACE_EXISTING);
					}

					//here we copy the input fasta file to the final output file
					else {
						Files.copy(FileSystems.getDefault().getPath(tmpFilePaths.get(0)),FileSystems.getDefault().getPath(outputPath),REPLACE_EXISTING);
					}
				}
				
				//paired end
				else {
					for(int i = 0; i < tmpFilePaths.size();i++) {
						if(hasNewIlluminaHeader(tmpFilePaths.get(i))) {
							tmpOutputFile = File.createTempFile("input", null, new File(tmpOutDir));
							renameNewIlluminaHeader(tmpFilePaths.get(i),tmpOutputFile.getAbsolutePath());
							tmpFilePaths.set(i, tmpOutputFile.getAbsolutePath());
						}
						
						else if(!hasPairedEndInfo(tmpFilePaths.get(i))) {
							tmpOutputFile = File.createTempFile("input", null, new File(tmpOutDir));
							addPairedEndInfo(tmpFilePaths.get(i),tmpOutputFile.getAbsolutePath(),(i+1));
							tmpFilePaths.set(i, tmpOutputFile.getAbsolutePath());
						}
					}
					
					concatenateFilesWithNIO(tmpFilePaths, outputPath);
				}
			}
			
			deleteFolderWithContent(new File(tmpOutDir));
		}
		catch(Exception e) {
			e.printStackTrace();
			deleteFolderWithContent(new File(tmpOutDir));
			System.exit(1);
		}
	}
	
	
	/**
	 * assumes two file paths in the list
	 * @param filePaths
	 * @return
	 */
	private boolean checkBaseNames(ArrayList<String> filePaths) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(new File(filePaths.get(0))));
		String firstHeader = br.readLine();
		br.close();
		
		br = new BufferedReader(new FileReader(new File(filePaths.get(1))));
		String secondHeader = br.readLine();
		br.close();
		return(firstHeader.substring(0,firstHeader.lastIndexOf('/')).equals(secondHeader.substring(0,firstHeader.lastIndexOf('/'))));
	}
	
	private FileType getFileType(String inputPath) {
		
		int i = inputPath.lastIndexOf('.');
		if(i != -1) {
			String extension = inputPath.substring(i + 1);
			if(extension.equals("fa") || extension.equals("fasta") || extension.equals("FA") || extension.equals("FASTA")) return FileType.FASTA;
			if(extension.equals("fastq") || extension.equals("Fastq") || extension.equals("FASTQ")) return FileType.FASTQ;
			if(extension.equals("zip") || extension.equals("Zip") || extension.equals("ZIP")) return FileType.ZIP;
			if(extension.equals("tar") || extension.equals("Tar") || extension.equals("TAR")) return FileType.TAR;
			if(extension.equals("gz") || extension.equals("Gz") || extension.equals("GZ")) {
				//check if it is tar.gz
				if(i > 2) {
					extension = inputPath.substring(i - 3);
					if(extension.equals("tar.gz") || extension.equals("TAR.GZ") || extension.equals("Tar.Gz")) return FileType.TARGZ;
					else return FileType.GZ;
				}
				else
					return FileType.GZ;
			}
			
			
			return FileType.UNKNOWN;
		}
		
		else {
			return FileType.UNKNOWN;
		}		
	}
	
	
	private ReadFormat getReadFormat(String readFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(readFilePath)));
			String header = br.readLine(); 
			if(header.charAt(0) == '>')
				return ReadFormat.FASTA;
			else if(header.charAt(0) == '@')
				return ReadFormat.FASTQ;
			else
				return ReadFormat.UNKNOWN;
		}
		catch(Exception e) {
			e.printStackTrace();
			return ReadFormat.UNKNOWN;
		}
	}
	
	private boolean hasNewIlluminaHeader(String fastaFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fastaFilePath)));
			String header = br.readLine().substring(1);
			br.close();
			if(header.contains(" ")) {
				String[] splittedHeader = header.split(" ")[1].split(":");
				if(splittedHeader.length >= 3 && (splittedHeader[0].equals("1") || splittedHeader[0].equals("2")) && (splittedHeader[1].equals("Y") || splittedHeader[1].equals("N"))) {
					return true;
				}
			}
			
			return false;
		}
		
		catch(Exception e) {
			e.printStackTrace();
			return false;
		}
	}
	
	private boolean hasPairedEndInfo(String fastaFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fastaFilePath)));
			String header = br.readLine().substring(1);
			br.close();
			if(header.length() > 2 && (header.substring(header.length() - 2).equals("/1") || header.substring(header.length() - 2).equals("/2"))) {
					return true;
			}
			
			return false;
		}
		
		catch(Exception e) {
			e.printStackTrace();
			return false;
		}
	}
	
	
	private void renameNewIlluminaHeader(String inputFilePath, String outputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedHeader;
			Pattern spacePattern = Pattern.compile(" ");
			Pattern doublePointPattern = Pattern.compile(":");
			
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					splittedHeader = spacePattern.split(currentLine);
					splittedHeader[0] += "/" + doublePointPattern.split(splittedHeader[1])[0];
					pw.println(splittedHeader[0]);
				}
				else {
					pw.println(currentLine);
				}
			}
			br.close();
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private void addPairedEndInfo(String inputFilePath, String outputFilePath, int mateInfo) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					currentLine += "/" + mateInfo;
				}
				
				pw.println(currentLine);
			}
			
			br.close();
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private void concatenateFilesWithNIO(ArrayList<String> filePaths, String outputPath) {
		try {
			File checkFile = new File(outputPath);
			if(checkFile.exists())
				checkFile.delete();
						
			FileOutputStream fos = new FileOutputStream(outputPath,true);
			FileChannel writeChannel = fos.getChannel();
			RandomAccessFile rf;
			FileChannel readChannel;
			long currentChannelSize;
			long transferedBytes;
			for(String filePath : filePaths) {
				if(!new File(filePath).exists())
					continue;
				
				rf = new RandomAccessFile(filePath,"r");
				readChannel = rf.getChannel();
				currentChannelSize = readChannel.size();
				transferedBytes = 0;
				while(transferedBytes < currentChannelSize) {
					transferedBytes += readChannel.transferTo(transferedBytes, readChannel.size(), writeChannel);
				}
				rf.close();
			}
			fos.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private void deleteFolderWithContent(File file) {
		if(file.isDirectory()) {
			File[] content = file.listFiles();
			for(File f : content)
				deleteFolderWithContent(f);
		}
		file.delete();
	}
	
	//TODO unzip all files in an archive (currently only the first file is unzipped)
	public void unzip(File inputFile, File outputFile) throws Exception {
		ZipInputStream inputStream = new ZipInputStream(new FileInputStream(inputFile));
		ZipEntry ze = inputStream.getNextEntry();
		
		//we expect that there is only a single file in the archive. If more fails are provided, only the first entry will be extracted...
		if(ze != null) {
			FileOutputStream outputStream = new FileOutputStream(outputFile);
			byte[] buffer = new byte[8192];
			int length;
			
			while((length=inputStream.read(buffer)) != -1) {
				outputStream.write(buffer, 0, length);
			}
			outputStream.close();
		}
		inputStream.closeEntry();
		inputStream.close();
		
	}
	
	
	public void gunzip(File inputFile, File outputFile) throws Exception {
		GZIPInputStream inputStream = new GZIPInputStream(new FileInputStream(inputFile));
		FileOutputStream outputStream = new FileOutputStream(outputFile);
		byte[] buffer = new byte[8192];
		int length;
		
		while((length=inputStream.read(buffer)) != -1) {
			outputStream.write(buffer, 0, length);
		}
		
		inputStream.close();
		outputStream.close();
	}
	
	
	
	
	public void convertFastqToFasta(File inputFile, File outputFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(inputFile));
			PrintWriter pw = new PrintWriter(new FileWriter(outputFile));
			
			String currentLine;
			String header;
			String sequence;
			
			while(br.ready()) {
				currentLine = br.readLine();
				
				if(currentLine.isEmpty())
					continue;
				
				if(currentLine.charAt(0) == '@') {
					header = String.format(">%s",currentLine.substring(1));
					sequence = br.readLine();
					pw.println(header);
					pw.println(sequence);
					//skip the next two fastq lines
					br.readLine();
					br.readLine();
				}
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
}
