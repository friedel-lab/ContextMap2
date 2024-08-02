package main;

import java.io.File;
import java.io.FileDescriptor;
import java.io.FileOutputStream;

import java.io.IOException;
import java.io.OutputStreamWriter;

public class UnsynchronizedFileWriter extends OutputStreamWriter {

	    public UnsynchronizedFileWriter(File file) throws IOException {
	        super(new UnsynchronizedFileOutputStream(file));
	    }

 
	
}
