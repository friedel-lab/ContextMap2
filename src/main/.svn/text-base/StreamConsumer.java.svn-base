package main;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class StreamConsumer extends Thread {
	
	private InputStream inputStream;
	private String outputType;
	
	private boolean verbose;
	private boolean processRunning;

	public StreamConsumer(InputStream inputStream, String outputType, boolean verbose) {
	this.inputStream = inputStream;
	this.outputType = outputType;
	this.verbose = verbose;
	this.processRunning = true;
	}
	
	public void processTerminated() {
		this.processRunning = false;
	}
	
	@Override
	public void run() {
		try {
			InputStreamReader inputStreamReader = new InputStreamReader(this.inputStream);
			BufferedReader br = new BufferedReader(inputStreamReader);
			while(!br.ready() && this.processRunning) {
				this.sleep(100);
			}
			
			//in case the stream is now ready, we output its content now
			if(br.ready()) {
				String line = null;
				while ((line = br.readLine()) != null) {
					if(this.outputType.equals("err")) {
						System.err.println(line);
					}
					else if(this.verbose) {
						System.out.println(line);
					}
				}
			}
			br.close();
			synchronized(this) {
				this.notifyAll();
			}
		}
		
		catch (Exception ioe) {
			ioe.printStackTrace();
			synchronized(this) {
				this.notifyAll();
			}
		}
	}
}
