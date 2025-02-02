/*
 * Copyright (c) 1996, 2005, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.  Oracle designates this
 * particular file as subject to the "Classpath" exception as provided
 * by Oracle in the LICENSE file that accompanied this code.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 */



package main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

public class UnsynchronizedBufferedWriter extends Writer {
	
	  private Writer out;

    private char cb[];
    private int nChars, nextChar;
    private static int defaultCharBufferSize = 8192;
    private String lineSeparator;

	public UnsynchronizedBufferedWriter(Writer out) {
		this(out,defaultCharBufferSize);
	}

	public UnsynchronizedBufferedWriter(Writer out, int sz) {
		super(out);
		 if (sz <= 0)
	            throw new IllegalArgumentException("Buffer size <= 0");
	        this.out = out;
	        cb = new char[sz];
	        nChars = sz;
	        nextChar = 0;

	        lineSeparator = System.getProperty("line.separator");
	}
	
	  /** Checks to make sure that the stream has not been closed */
    private void ensureOpen() throws IOException {
        if (this.out == null)
            throw new IOException("Stream closed");
    }


    /**
     * Flushes the output buffer to the underlying character stream, without
     * flushing the stream itself.  This method is non-private only so that it
     * may be invoked by PrintStream.
     */
    public void flushBuffer() throws IOException {
        ensureOpen();
        if (nextChar == 0)
            return;
        out.write(cb, 0, nextChar);
        nextChar = 0;
    }

    /**
     * Writes a single character.
     *
     * @exception  IOException  If an I/O error occurs
     */
    public void write(int c) throws IOException {
        ensureOpen();
        if (nextChar >= nChars)
            flushBuffer();
        cb[nextChar++] = (char) c;
    }

    /**
     * Our own little min method, to avoid loading java.lang.Math if we've run
     * out of file descriptors and we're trying to print a stack trace.
     */
    private int min(int a, int b) {
        if (a < b) return a;
        return b;
    }

    /**
     * Writes a portion of an array of characters.
     *
     * <p> Ordinarily this method stores characters from the given array into
     * this stream's buffer, flushing the buffer to the underlying stream as
     * needed.  If the requested length is at least as large as the buffer,
     * however, then this method will flush the buffer and write the characters
     * directly to the underlying stream.  Thus redundant
     * <code>BufferedWriter</code>s will not copy data unnecessarily.
     *
     * @param  cbuf  A character array
     * @param  off   Offset from which to start reading characters
     * @param  len   Number of characters to write
     *
     * @exception  IOException  If an I/O error occurs
     */
    public void write(char cbuf[], int off, int len) throws IOException {
        ensureOpen();
        if ((off < 0) || (off > cbuf.length) || (len < 0) ||
            ((off + len) > cbuf.length) || ((off + len) < 0)) {
            throw new IndexOutOfBoundsException();
        } else if (len == 0) {
            return;
        }

        if (len >= nChars) {
            /* If the request length exceeds the size of the output buffer,
               flush the buffer and then write the data directly.  In this
               way buffered streams will cascade harmlessly. */
            flushBuffer();
            out.write(cbuf, off, len);
            return;
        }

        int b = off, t = off + len;
        while (b < t) {
            int d = min(nChars - nextChar, t - b);
            System.arraycopy(cbuf, b, cb, nextChar, d);
            b += d;
            nextChar += d;
            if (nextChar >= nChars)
                flushBuffer();
        }
    }

    /**
     * Writes a portion of a String.
     *
     * <p> If the value of the <tt>len</tt> parameter is negative then no
     * characters are written.  This is contrary to the specification of this
     * method in the {@linkplain java.io.Writer#write(java.lang.String,int,int)
     * superclass}, which requires that an {@link IndexOutOfBoundsException} be
     * thrown.
     *
     * @param  s     String to be written
     * @param  off   Offset from which to start reading characters
     * @param  len   Number of characters to be written
     *
     * @exception  IOException  If an I/O error occurs
     */
    public void write(String s, int off, int len) throws IOException {
        ensureOpen();

        int b = off, t = off + len;
        while (b < t) {
            int d = min(nChars - nextChar, t - b);
            s.getChars(b, b + d, cb, nextChar);
            b += d;
            nextChar += d;
            if (nextChar >= nChars)
                flushBuffer();
        }
    }

    /**
     * Writes a line separator.  The line separator string is defined by the
     * system property <tt>line.separator</tt>, and is not necessarily a single
     * newline ('\n') character.
     *
     * @exception  IOException  If an I/O error occurs
     */
    public void newLine() throws IOException {
        write(lineSeparator);
    }

    /**
     * Flushes the stream.
     *
     * @exception  IOException  If an I/O error occurs
     */
    
    public void flush() throws IOException {
        flushBuffer();
        out.flush();
    }

    public void close() throws IOException {
            if (out == null) {
                return;
            }
            try {
                flushBuffer();
            } finally {
                out.close();
                out = null;
                cb = null;
            }
    }

}
