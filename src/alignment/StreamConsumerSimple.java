package alignment;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;


public class StreamConsumerSimple extends SamStreamConsumer {
	
	private final String OUTPUT_FILE_PATH;	
	private static final String NEWLINE = System.lineSeparator();

	public StreamConsumerSimple(InputStream inputStream, StreamType streamType, String outputFileName) {
		super(inputStream, streamType);
		this.OUTPUT_FILE_PATH = outputFileName;
	}
	
	@Override
	public void run() {
		try {
			// open the reader and writer
			BufferedWriter bw = null;
			if(this.OUTPUT_FILE_PATH != null)
				bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.OUTPUT_FILE_PATH), true)));
				
			InputStreamReader inputStreamReader = new InputStreamReader(this.inputStream);
			BufferedReader br = new BufferedReader(inputStreamReader);
			// wait until buffer is full or thread is not running any more
			while(!br.ready() && this.processRunning.get()) {
				sleep(100);
			}
			
			if(br.ready()) {
				String line;
				// read the lines and write them, if needed
				while((line = br.readLine()) != null) {
					if(bw != null) {
						bw.write(line);
						bw.append(NEWLINE);
					}
				}
			}
			// close the buffer and flush 
			br.close();
			if(bw != null) {
				bw.flush();
				bw.close();
			}
			
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