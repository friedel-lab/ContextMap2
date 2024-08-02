package alignment;

import java.io.InputStream;
import java.util.concurrent.atomic.AtomicBoolean;

public class SamStreamConsumer extends Thread {
	
	protected InputStream inputStream;
	protected StreamType streamType;
	protected AtomicBoolean processRunning;
	
	public SamStreamConsumer(InputStream inputStream, StreamType streamType) {
		this.inputStream = inputStream;
		this.streamType  = streamType;
		this.processRunning = new AtomicBoolean(true);
	}
	
	public void processTerminated() {
		this.processRunning.set(false);
	}
}
