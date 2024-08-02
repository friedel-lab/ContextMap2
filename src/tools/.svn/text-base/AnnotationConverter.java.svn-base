package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

public class AnnotationConverter {

	
	public static void main(String args[]) {
		if(args[0].equals("gtf2cm")) {
			String inputFilePath = null;
			String customChrName = null;
			boolean printCDS = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					inputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-customChrName")) {
					customChrName = args[++i];
					continue;
				}
				
				if(args[i].equals("--cds")) {
					printCDS = true;
					continue;
				}
				
			}
			gtf2cm(inputFilePath,customChrName,printCDS);
		}
		
		if(args[0].equals("csv2cm")) {
			String inputFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					inputFilePath = args[++i];
					continue;
				}
				
			}
			csv2cm(inputFilePath);
		}
		
		if(args[0].equals("remapHsvIds")) {
			String gtfFilePath = null;
			String cmAnnotationFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-gtf")) {
					gtfFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-cm")) {
					cmAnnotationFilePath = args[++i];
					continue;
				}
				
			}
			remapHsvGeneNames(gtfFilePath,cmAnnotationFilePath);
		}
		
		if(args[0].equals("cm2geneList")) {
			String outputFilePath = null;
			String cmAnnotationFilePath = null;
			String biotype = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-i")) {
					cmAnnotationFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-biotype")) {
					biotype = args[++i];
					continue;
				}
				
			}
			cm2geneList(cmAnnotationFilePath,outputFilePath,biotype);
		}
	}
	
	
	private static void cm2geneList(String cmFilePath, String outputFilePath, String biotype) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(cmFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			
			String currentLine;
			String[] splittedLine;
			Pattern spacePattern = Pattern.compile(" ");
			pw.println("gene_id,chr_name,biotype");
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '#') {
					splittedLine = spacePattern.split(currentLine);
					pw.println(String.format("%s,%s,%s",splittedLine[1],splittedLine[2],biotype));
				}
			}
			br.close();
			pw.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void remapHsvGeneNames(String gtfFilePath, String cmAnnotationFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(gtfFilePath)));
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern semiColonPattern = Pattern.compile(";");
			String xref;
			String id;
			HashMap<String,String> xref2id = new HashMap<String,String>();
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				if(splittedLine[8].contains("db_xref")) {
					splittedLine = semiColonPattern.split(splittedLine[8]);
					id = splittedLine[0].substring(splittedLine[0].indexOf('"') + 1,splittedLine[0].lastIndexOf('"'));
					xref = splittedLine[3].substring(splittedLine[3].indexOf('"') + 1,splittedLine[3].lastIndexOf('"'));
					if(!xref2id.containsKey(id)) {
						xref2id.put(id, xref);}
					else {
						if(!xref2id.get(id).equals(xref)) {
							System.out.println("Found more than one possible mapping for gene id: " + id);
							System.exit(1);
						}
					}
					
				}
			}
			br.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void csv2cm(String inputFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));
			//skip header
			br.readLine();
			
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			
			String geneId;
			String prevGeneId = "";
			Gene currentGene = null;
			
			String transcriptId;
			String prevTranscriptId = "";
			Transcript transcript = null;
			
			String exonId;
			String chr;
			String strand;
			int start;
			int end;
			Exon exon = null;
			
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				geneId = splittedLine[0];
				transcriptId = splittedLine[1];
				exonId = splittedLine[2];
				chr = splittedLine[3];
				strand = splittedLine[4].equals("1")?"+":"-";
				start = Integer.valueOf(splittedLine[5]);
				end = Integer.valueOf(splittedLine[6]);
				exon = new Exon(exonId,start,end);
				
				if(!geneId.equals(prevGeneId)) {
					if(currentGene != null) {
						writeGene(currentGene);
					}
					currentGene = new Gene(geneId,chr,strand,start,end);
				}
				
					
				if(!transcriptId.equals(prevTranscriptId)) {
					transcript = new Transcript(transcriptId,start,end);
					currentGene.addTranscript(transcript);
				}
				
				if(!currentGene.containsExon(exon.getId())) {
					currentGene.addExon(exon);
					if(currentGene.getEnd() < end)
						currentGene.setEnd(end);
					if(currentGene.getStart() > start)
						currentGene.setStart(start);
					
				}
				
				transcript.addExon(currentGene.getExonPosition(exon.getId()));
				if(transcript.getEnd() < end)
					transcript.setEnd(end);
				if(transcript.getStart() > start) 
					transcript.setStart(start);

				
				
				prevGeneId = geneId;
				prevTranscriptId = transcriptId;
			}
			
			br.close();
			
			if(currentGene != null) {
				writeGene(currentGene);
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}

	
	private static void gtf2cm(String inputFilePath, String customChrName, boolean printCDS) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(inputFilePath)));			
			String chr;
			String strand;
			Gene currentGene = null;
			String currentGeneId = null;
			String prevGeneId = "";
			
			Transcript currentTranscript = null;
			String currentTranscriptId = null;
			String prevTranscriptId = "";
			
			
			Exon currentExon;
			String exonId = null;
			int start;
			int end;
			String extendedInfo;
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern semicolonPattern = Pattern.compile(";");
			
			while((currentLine = br.readLine()) != null) {
				
				if(currentLine.charAt(0) == '#')
					continue;
				
				splittedLine = tabPattern.split(currentLine);
				if(!splittedLine[2].equals("exon") && !(printCDS && splittedLine[2].equals("CDS")))
					continue;
				
				chr = splittedLine[0];
				if(chr.length() < 3 || !chr.substring(0,3).equals("chr"))
					chr = "chr" + chr;
				
				if(customChrName != null) {
					chr = customChrName;
				}
				
				start = Integer.valueOf(splittedLine[3]);
				end = Integer.valueOf(splittedLine[4]);
				strand = splittedLine[6];
				extendedInfo = splittedLine[8];
				exonId = null;
				splittedLine = semicolonPattern.split(extendedInfo);
				for(int i = 0; i < splittedLine.length; i++) {
					if(splittedLine[i].contains("gene_id")) {
						currentGeneId = splittedLine[i].substring(splittedLine[i].indexOf('"') + 1,splittedLine[i].lastIndexOf('"'));
						continue;
					}
					
					if(splittedLine[i].contains("transcript_id")) {
						currentTranscriptId = splittedLine[i].substring(splittedLine[i].indexOf('"') + 1,splittedLine[i].lastIndexOf('"'));
						continue;
					}
					
					if(splittedLine[i].contains("exon_id")) {
						exonId = splittedLine[i].substring(splittedLine[i].indexOf('"') + 1,splittedLine[i].lastIndexOf('"'));
					}
				}
				
				if(!currentGeneId.equals(prevGeneId)) {
					if(currentGene != null) {
						writeGene(currentGene);
					}
					
					currentGene = new Gene(currentGeneId,chr,strand,start,end);
				}
				
				if(!currentTranscriptId.equals(prevTranscriptId)) {
					currentTranscript = new Transcript(currentTranscriptId,start,end);
					currentGene.addTranscript(currentTranscript);
				}
				
				if(exonId == null) {
					exonId = String.format("%s_%s",start,end);
				}
				
				
				if(!currentGene.containsExon(exonId)) {
					currentExon = new Exon(exonId,start,end);
					currentGene.addExon(currentExon);
					if(currentGene.getEnd() < end)
						currentGene.setEnd(end);
					if(currentGene.getStart() > start)
						currentGene.setStart(start);
					
				}
				currentTranscript.addExon(currentGene.getExonPosition(exonId));
				if(currentTranscript.getEnd() < end)
					currentTranscript.setEnd(end);
				if(currentTranscript.getStart() > start) 
					currentTranscript.setStart(start);
				
				
				
				prevGeneId = currentGeneId;
				prevTranscriptId = currentTranscriptId;
			}
			br.close();
			
			if(currentGene != null) {
				writeGene(currentGene);
			}
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void writeGene(Gene g) {
		System.out.println(String.format("# %s %s %s %s %s",g.getId(),g.getChr(),g.getStrand(),g.getStart(),g.getEnd()));
		ArrayList<Exon> exons = g.getExons();
		System.out.print(String.format("%s %s %s",exons.get(0).getId(),exons.get(0).getStart(),exons.get(0).getEnd()));
		for(int i = 1; i < exons.size(); i++) {
			System.out.print(String.format("\t%s %s %s",exons.get(i).getId(),exons.get(i).getStart(),exons.get(i).getEnd()));
		}
		System.out.println();
		ArrayList<Transcript> transcripts = g.getTranscripts();
		ArrayList<Integer> exonNumbers = transcripts.get(0).getExons();
		System.out.print(String.format("%s %s %s %s",transcripts.get(0).getId(),transcripts.get(0).getStart(),transcripts.get(0).getEnd(),exonNumbers.get(0)));
		for(int j = 1; j < exonNumbers.size(); j++) {
			System.out.print("," + exonNumbers.get(j));
		}
		
		for(int i = 1; i < transcripts.size(); i++) {
			exonNumbers = transcripts.get(i).getExons();
			System.out.print(String.format("\t%s %s %s %s",transcripts.get(i).getId(),transcripts.get(i).getStart(),transcripts.get(i).getEnd(),exonNumbers.get(0)));
			for(int j = 1; j < exonNumbers.size(); j++) {
				System.out.print("," + exonNumbers.get(j));
			}
		}
		System.out.println();
	}
	
	
	
	private static class Gene {
		private String id;
		private String chr;
		private String strand;
		private int start;
		private int end;
		
		private HashMap<String, Integer> exonId2arrayPosition;
		private ArrayList<Exon> exons;
		private ArrayList<Transcript> transcripts;
		
		public Gene(String id, String chr, String strand, int start, int end) {
			this.id = id;
			this.chr = chr;
			this.strand = strand;
			this.start = start;
			this.end = end;
			
			this.exons = new ArrayList<Exon>();
			this.exonId2arrayPosition = new HashMap<String,Integer>();
			this.transcripts = new ArrayList<Transcript>();
		}

		public String getId() {
			return id;
		}

		public String getChr() {
			return chr;
		}

		public String getStrand() {
			return strand;
		}

		public int getStart() {
			return start;
		}
		
		public void setStart(int start) {
			this.start  = start;
		}

		public int getEnd() {
			return end;
		}
		
		public void setEnd(int end) {
			this.end = end;
		}
		
		public ArrayList<Exon> getExons() {
			return this.exons;
		}
		
		public void addExon(Exon e) {
			this.exons.add(e);
			this.exonId2arrayPosition.put(e.getId(),this.exons.size() - 1);
		}
		
		public boolean containsExon(String id) {
			return this.exonId2arrayPosition.containsKey(id);
		}
		
		public int getExonPosition(String id) {
			return this.exonId2arrayPosition.get(id);
		}
		
		public ArrayList<Transcript> getTranscripts() {
			return this.transcripts;
		}
		
		public void addTranscript(Transcript t) {
			this.transcripts.add(t);
		}
	}
	
	
	private static class Exon {
		private String id;
		private int start;
		private int end;
		
		public Exon(String id, int start, int end) {
			this.id = id;
			this.start = start;
			this.end = end;
		}
		
		public String getId() {
			return id;
		}
		
		public int getStart() {
			return start;
		}

		public int getEnd() {
			return end;
		}
	}
	
	private static class Transcript {
		private String id;
		private int start;
		private int end;
		private ArrayList<Integer> exons;
		
		public Transcript(String id, int start, int end) {
			this.id = id;
			this.start = start;
			this.end = end;
			this.exons = new ArrayList<Integer>();
		}
		
		public String getId() {
			return id;
		}
		
		public int getStart() {
			return start;
		}
		
		public void setStart(int start) {
			this.start = start;
		}

		public int getEnd() {
			return end;
		}
		
		public void setEnd(int end) {
			this.end = end;
		}
		
		public void addExon(int i) {
			this.exons.add(i);
		}
		
		public ArrayList<Integer> getExons() {
			return this.exons;
		}
		
	} 
	
}
