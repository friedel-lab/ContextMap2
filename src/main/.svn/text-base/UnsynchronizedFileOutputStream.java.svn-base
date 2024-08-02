package main;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;


public class UnsynchronizedFileOutputStream extends FileOutputStream {

	 private volatile boolean closed = false;
	
	public UnsynchronizedFileOutputStream(File file) throws FileNotFoundException {
		super(file);
	}
   
	public void close() throws IOException {
            if (this.closed) {
                return;
            }
            closed = true;
    }
	
}
