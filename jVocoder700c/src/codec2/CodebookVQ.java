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
import java.util.ArrayList;
import java.util.Scanner;

public final class CodebookVQ implements IDefines {

    private Scanner m_scanner;
    private InputStreamReader m_stream;
    private final ArrayList<Float> m_codes0;
    private final ArrayList<Float> m_codes1;
    
    protected record Codebook(ArrayList<Float> CodeBookArray) {

        public float CodeBookArraySubscript(int index) {
            return CodeBookArray.get(index);
        }
    }

    private final Codebook[] m_vqcb;

    /*
     * We have to read these VQ values in, as the number of values is too
     * great for the Java compiler to handle as static constants.
     *
     * These text files are set as resources in the JAR file, so the user
     * doesn't have to maintain them. They are in the package 'codebook'
     */
    protected CodebookVQ() {
        m_vqcb = new Codebook[2];

        try {
            m_stream = new InputStreamReader(getClass().getResourceAsStream("/codebook/train_120_1.txt"));
            m_scanner = new Scanner(new BufferedReader(m_stream));
        } catch (Exception fnfe1) {
            System.err.println("Vector file not found '/codebook/train_120_1.txt' " + fnfe1.getMessage());
            System.exit(1);
        }

        m_codes0 = new ArrayList<>(AMP_K * AMP_M); // 10,240 (k * m)

        int index = 0;
        while (m_scanner.hasNextFloat()) {
            m_codes0.add(index, m_scanner.nextFloat());
            index++;
        }

        m_vqcb[0] = new Codebook(m_codes0); // k, log2, m

        m_scanner.close();

        try {
            m_stream = new InputStreamReader(getClass().getResourceAsStream("/codebook/train_120_2.txt"));
            m_scanner = new Scanner(new BufferedReader(m_stream));
        } catch (Exception fnfe2) {
            System.err.println("Vector file not found '/codebook/train_120_2.txt' " + fnfe2.getMessage());
            System.exit(2);
        }

        m_codes1 = new ArrayList<>(AMP_K * AMP_M); // 10,240 (k * m)

        index = 0;
        while (m_scanner.hasNextFloat()) {
            m_codes1.add(index, m_scanner.nextFloat());
            index++;
        }
        
        m_vqcb[1] = new Codebook(m_codes1); // k, log2, m
        
        m_scanner.close();
    }

    /**
     * Method used to reference codebook by amp.java
     *
     * @param index int which is an index reference
     * @return codebook vector quantized codebook
     */
    protected Codebook getCodebook(int index) {
        return m_vqcb[index];
    }
}
