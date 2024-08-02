package tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;

public class String2Bitset {

	/**
	 * A -> 001
	 * C -> 010
	 * G -> 011
	 * T -> 100
	 * N -> 101
	 */
	
	private HashMap<Character,char[]> base2bits = new HashMap<Character,char[]>();
	private HashMap<String,Character> bits2base = new HashMap<String,Character>();
	
	private long[] prevArray;
	private String prevSequence;
	
	public String2Bitset() {
	
		base2bits = new HashMap<Character,char[]>();
		base2bits.put('A',new char[]{'0','0','1'});
		base2bits.put('C',new char[]{'0','1','0'});
		base2bits.put('G',new char[]{'0','1','1'});
		base2bits.put('T',new char[]{'1','0','0'});
		base2bits.put('N',new char[]{'1','0','1'});
		
		bits2base.put("001", 'A');
		bits2base.put("010", 'C');
		bits2base.put("011", 'G');
		bits2base.put("100", 'T');
		bits2base.put("101", 'N');
		
		this.prevArray = null;
		this.prevSequence = null;
	}
	
	
	 public long[] compress(String s) {
	        BitSet bs = new BitSet(s.length() * 3);
	        for (int i = 0, j = 0; i < s.length(); i++, j+=3) {
	        	char[] v = this.base2bits.get(s.charAt(i));
	            for (int k = j, l = 0; l < 3; l++, k++) {
	                if ('0' == v[l])
	                    continue;
	                
	                bs.set(k);
	            }
	        }
	        return bs.toLongArray();
	 }
	 
	 
	 public String decompress(long[] array) {
		 StringBuilder tmpBuffer = new StringBuilder(3);
		 BitSet bs = BitSet.valueOf(array);
		 StringBuilder sb = new StringBuilder(bs.length()/3);
		 for(int i = 0; i < bs.length(); i+= 3) {
			 tmpBuffer.setLength(0);
			 for(int j = i; j < i+3; j++) {
				 if(bs.get(j))
					 tmpBuffer.append('1');
				 else
					 tmpBuffer.append('0');
				 }
			 sb.append(this.bits2base.get(tmpBuffer.toString()));
		 }
		 
		 return sb.toString();
	 }
	 
	 
/*	 public String decompressBuffered(long[] array) {
		 
		 if(Arrays.equals(array, this.prevArray))
			 return this.prevSequence;
		 
		 StringBuilder tmpBuffer = new StringBuilder(3);
		 BitSet bs = BitSet.valueOf(array);
		 StringBuilder sb = new StringBuilder(bs.length()/3);
		 for(int i = 0; i < bs.length(); i+= 3) {
			 tmpBuffer.setLength(0);
			 for(int j = i; j < i+3; j++) {
				 if(bs.get(j))
					 tmpBuffer.append('1');
				 else
					 tmpBuffer.append('0');
				 }
			 sb.append(this.bits2base.get(tmpBuffer.toString()));
		 }
		 
		 this.prevArray = array;
		 this.prevSequence = sb.toString();
		 return this.prevSequence;
	 }
*/	
}

