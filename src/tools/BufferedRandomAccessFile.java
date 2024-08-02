package tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;

public class BufferedRandomAccessFile extends RandomAccessFile {

	private byte buffer[];
	private int bufferSize;
	private int bufEnd;
	private int bufPos;
	private long realPos;
	
	public BufferedRandomAccessFile(File file, String mode, int bufferSize) throws IOException {
		super(file, mode);
		invalidate();
		this.bufferSize = bufferSize;
		this.buffer = new byte[bufferSize];
	}
	
	public final int read() throws IOException {
		if(this.bufPos >= this.bufEnd) {
			if(fillBuffer() < 0)
				return -1;
		}
		
		if(this.bufEnd == 0)
			return -1;
		 
		else
			return this.buffer[this.bufPos++];
	}

	private int fillBuffer() throws IOException {
	    int n = super.read(this.buffer, 0, this.bufferSize);
	    if(n >= 0) {
	      this.realPos +=n;
	      this.bufEnd = n;
	      this.bufPos = 0;
	    }
	    return n;
	}
	
	private void invalidate() throws IOException {
		this.bufEnd = 0;
		this.bufPos = 0;
		this.realPos = super.getFilePointer();
	  }
	
	public int read(byte b[], int off, int length) throws IOException {
		int leftover = this.bufEnd - this.bufPos;
		if(length <= leftover) {
			System.arraycopy(this.buffer, this.bufPos, b, off, length);
			this.bufPos += length;
			return length;
		}
		for(int i = 0; i < length; i++) {
			int c = this.read();
			if(c != -1)
				b[off+i] = (byte)c;
			else {
				return i;
			}
		}
		return length;
	}
	
	public long getFilePointer() throws IOException{
	    long l = this.realPos;
	    return (l - this.bufEnd + this.bufPos) ;
	  }
	
	  public void seek(long pos) throws IOException {
	    long n = this.realPos - pos;
	    if(n >= 0 && n <= this.bufEnd) {
	    	this.bufPos = this.bufEnd - (int)n;
	    } 
	    else {
	      super.seek(pos);
	      invalidate();
	    }
	  }
	  
	  public final String getNextLine() throws IOException {
		  String str = null;
		  if(this.bufEnd - this.bufPos <= 0) {
			  if(fillBuffer() < 0) {
				  return null;
			  }
		  }
		  int lineEnd = -1;
		  for(int i = this.bufPos; i < this.bufEnd; i++) {
			  if(buffer[i] == '\n') {
				  lineEnd = i;
				  break;
			  }
		  }
	   if(lineEnd < 0) {
	        StringBuilder input = new StringBuilder(256);
	        int c;
	             while (((c = read()) != -1) && (c != '\n')) {
	                 input.append((char)c);
	        }
	        if ((c == -1) && (input.length() == 0)) {
	          return null;
	        }
	        return input.toString();
	   }
	   if(lineEnd > 0 && this.buffer[lineEnd-1] == '\r')
	        str = new String(this.buffer, 0, this.bufPos, lineEnd - this.bufPos -1);
	   else str = new String(this.buffer, 0, this.bufPos, lineEnd - this.bufPos);
	   this.bufPos = lineEnd + 1;
	   return str;
	  }
	 
	
}
