/*
 * Copyright (C) 1993-2017 David Rowe
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 *
 * Java Version by Steve Sampson
 */
package codec2;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Scanner;

/*
 * I'm using floats here, to save memory space (4 versus 8 bytes each)
 */
public final class CodebookVQ {

    private Scanner m_scanner;
    private InputStreamReader m_stream;
    private final Codebook[] m_vqcb;
    private final float[] m_codes0;
    private final float[] m_codes1;

    /*
     * We have to read these VQ values in, as the number of values is too
     * great for the Java compiler to handle as static constants.
     *
     * These text files are set as resources in the JAR file, so the user
     * doesn't have to maintain them. They are in the package 'codebook'
     */
    public CodebookVQ() {
        m_vqcb = new Codebook[2];

        try {
            m_stream = new InputStreamReader(getClass().getResourceAsStream("/codebook/train_120_1.txt"));
            m_scanner = new Scanner(new BufferedReader(m_stream));
        } catch (Exception fnfe1) {
            System.err.println("Vector file not found '/codebook/train_120_1.txt' " + fnfe1.getMessage());
            System.exit(1);
        }

        int k = m_scanner.nextInt();
        int m = m_scanner.nextInt();
        int log2m = (int) (Math.log((double) m) / Math.log(2.0));

        m_codes0 = new float[k * m];
        m_vqcb[0] = new Codebook(k, log2m, m, m_codes0);

        int index = 0;
        while (m_scanner.hasNextFloat()) {
            m_codes0[index] = m_scanner.nextFloat();
            index++;
        }

        m_scanner.close();

        try {
            m_stream = new InputStreamReader(getClass().getResourceAsStream("/codebook/train_120_2.txt"));
            m_scanner = new Scanner(new BufferedReader(m_stream));
        } catch (Exception fnfe2) {
            System.err.println("Vector file not found '/codebook/train_120_2.txt' " + fnfe2.getMessage());
            System.exit(2);
        }

        k = m_scanner.nextInt();
        m = m_scanner.nextInt();
        log2m = (int) (Math.log((double) m) / Math.log(2.0));

        m_codes1 = new float[k * m];
        m_vqcb[1] = new Codebook(k, log2m, m, m_codes1);

        index = 0;
        while (m_scanner.hasNextFloat()) {
            m_codes1[index] = m_scanner.nextFloat();
            index++;
        }

        m_scanner.close();
    }

    /**
     * Method used to reference codebook by amp.java
     *
     * @param index int which is an index reference
     * @return codebook vector quantized codebook
     */
    public Codebook getCodebook(int index) {
        return m_vqcb[index];
    }
}
