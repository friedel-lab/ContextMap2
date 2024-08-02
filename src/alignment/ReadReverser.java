package alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

public class ReadReverser {

	public ReadReverser() {
		
	}
	
	
	public void getReverseComplement(String inputFilePath, String outputFilePath, String readFormat, boolean addRcInfo) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			String currentLine;
			StringBuilder sb = new StringBuilder();
			
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>' || currentLine.charAt(0) == '@') {
					pw.print(currentLine);
					if(addRcInfo)
						pw.print("/rc");
					pw.println();
					
					sb.setLength(0);
					sb.append(br.readLine());
					sb.reverse();
					for(int i = 0; i < sb.length(); i++) {
						sb.setCharAt(i, substitute(sb.charAt(i)));
					}
					pw.println(sb.toString());
					
					if(readFormat.equals("fastq")) {
						pw.println(br.readLine());
						sb.setLength(0);
						sb.append(br.readLine());
						pw.println(sb.reverse().toString());
					}
				}
			}
			br.close();
			pw.close();
		}
		
		catch(Exception e ) {
			e.printStackTrace();
		}
	}
	
	
	
	private char substitute(char n) {
		if(n == 'A' || n == 'a') return 'T';
		if(n == 'T' || n == 't') return 'A';
		if(n == 'C' || n == 'c') return 'G';
		if(n == 'G' || n == 'g') return 'C';
		else return 'N';
	}
}
