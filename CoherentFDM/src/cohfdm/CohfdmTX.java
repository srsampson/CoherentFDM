/*
 * Copyright (C) 2015 David Rowe
 *
 * All rights reserved
 *
 * Java Version by Steve Sampson
 */
package cohfdm;

import complex.Complex;
import complex.ComplexMath;

/**
 * Coherent FDM Modem Transmitter Class
 */
public final class CohfdmTX implements IDefines {

    private final Complex[][] txFilterMemory;
    private final Complex[] mtxPhase;
    private final Complex[] mtxFrequency;
    private Complex mbasebandTxPhase;
    private float mfrequencySeparation;
    private int mdiversityFactor;

    /**
     * Instantiate with Diversity True
     */
    public CohfdmTX() {
        this(true);
    }

    /**
     * Instantiate with Diversity option true = Diversity, false = No Diversity
     *
     * @param div boolean to enable or disable diversity
     */
    public CohfdmTX(boolean div) {
        int r, c;

        mdiversityFactor = (div) ? 2 : 1;

        mfrequencySeparation = 75.0f * 1.5f;    // 75 Hz + .5 Excess Bandwidth

        txFilterMemory = new Complex[COHPSK_NC * mdiversityFactor][COHPSK_NSYM];
        mtxPhase = new Complex[COHPSK_NC * mdiversityFactor];
        mtxFrequency = new Complex[COHPSK_NC * mdiversityFactor];

        for (c = 0; c < (COHPSK_NC * mdiversityFactor); c++) {
            for (r = 0; r < COHPSK_NSYM; r++) {
                txFilterMemory[c][r] = new Complex();
            }
        }

        Complex val = ComplexMath.cexp(new Complex());
        mbasebandTxPhase = val;

        for (c = 0; c < (COHPSK_NC * mdiversityFactor); c++) {
            mtxPhase[c] = val;

            /*
             * Spread the frequencies out by different spreading values
             */
            float freqHz = mfrequencySeparation * (-(COHPSK_NC * mdiversityFactor) / 2.0f - 0.5f + (float) Math.pow(c + 1, 0.98));

            mtxFrequency[c] = ComplexMath.cexp(new Complex(0.0f, TAU * freqHz / COHPSK_FS));
        }
    }

    public boolean getDiversity() {
        return (mdiversityFactor == 2);
    }

    public float getFreqSeparation() {
        return mfrequencySeparation;
    }

    /**
     * Method to return a modulated signal of the given bits
     *
     * @param symbols a complex array of the modulated signal
     * @param bits a boolean array of the data bits
     */
    public void modulate(Complex[] symbols, boolean[] bits) {
        Complex[][] tx_symb = new Complex[NSYMROWPILOT][COHPSK_NC * mdiversityFactor];
        Complex[] tx_onesym = new Complex[COHPSK_NC * mdiversityFactor];

        bitsToSymbols(tx_symb, bits);

        // create the 6 row by QPSK column modem frame
        for (int r = 0; r < NSYMROWPILOT; r++) {
            System.arraycopy(tx_symb[r], 0, tx_onesym, 0, COHPSK_NC * mdiversityFactor);

            upconvert(symbols, (r * COHPSK_M), tx_onesym);
        }

        // Hilbert Clipping to reduce Crest Factor by about 2 dB
        // this will typically occur about 5% of the signal samples
        for (int r = 0; r < COHPSK_NOM_TX_SAMPLES_PER_FRAME; r++) {
            float mag = 1E-12f + ComplexMath.cabs(symbols[r]); // divide by zero protect

            if (mag > COHPSK_CLIP) {        // 6.5
                symbols[r] = ComplexMath.times(symbols[r], COHPSK_CLIP / mag);
            }
        }
    }

    /**
     * QPSK Modulate the data and pilot bits
     *
     * @param symbols the modulated pilot and data bits
     * @param bits the data bits to be modulated
     */
    private void bitsToSymbols(Complex[][] symbols, boolean[] bits) {
        int row, column, newcolumn, pilotRow, dataRow, bitPair;

        /*
         * Organize QPSK symbols into 6 rows by 7 QPSK column frames.
         * Then duplicate these to create 14 QPSK carriers total.
         *
         * Insert 2 rows of pilots at beginning of frame. (28 bits)
         * Insert 4 rows of data after the pilots.        (56 bits)
         *
         * Each column is a carrier, each of the 6 rows (84 bits total)
         * is in time order (13.33 ms per row).
         *
         * sqrt(ND) term ensures the same energy/symbol for different
         * diversity factors.
         */
        for (row = 0, pilotRow = 0; pilotRow < 2; pilotRow++, row++) {
            for (column = 0; column < COHPSK_NC; column++) {
                if (mdiversityFactor == 1) {      // save some CPU from a complex divide
                    symbols[row][column] = new Complex(PILOTS[pilotRow][column], 0.0f);
                } else {
                    symbols[row][column] = new Complex(PILOTS[pilotRow][column] / 1.414213562f, 0.0f); // sqrt 2.0
                }
            }
        }

        for (dataRow = 0; dataRow < NSYMROW; dataRow++, row++) {
            for (column = 0; column < COHPSK_NC; column++) {
                int i = column * NSYMROW + dataRow;
                bitPair = ((bits[2 * i] ? 1 : 0) << 1) | (bits[2 * i + 1] ? 1 : 0);

                if (mdiversityFactor == 1) {      // save some CPU from a complex divide
                    symbols[row][column] = CONSTELLATION[bitPair];
                } else {
                    symbols[row][column] = ComplexMath.divide(CONSTELLATION[bitPair], 1.414213562f); // sqrt 2.0
                }
            }
        }

        /*
         * Duplicate the 7 QPSK carriers for diversity if enabled
         */
        if (mdiversityFactor == 2) {
            for (row = 0; row < NSYMROWPILOT; row++) {
                for (column = 0, newcolumn = COHPSK_NC; newcolumn < (COHPSK_NC * 2); column++, newcolumn++) {
                    symbols[row][newcolumn] = symbols[row][column];
                }
            }
        }
    }

    /**
     * Mix the new baseband signal segment into the output waveform.
     *
     * The output waveform is 600 symbols, and each baseband segment is 100
     * symbols, so the offset variable is used to position the mixed signal in
     * the proper position of the waveform array.
     *
     * This is basically a Java thing, as you can't use pointers as in the C
     * macro assembler language types.
     *
     * @param waveform the complex output signal centered on 1500 Hz
     * @param offset the index to put the mixed signal in the waveform array
     * @param baseband the input baseband signal to be mixed up
     */
    private void upconvert(Complex[] waveform, int offset, Complex[] baseband) {
        /*
         * Add the (baseband symbols * .707 amplitude) to the last row,
         * which is available now since the earlier call has
         * moved each row up one row.
         */
        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            txFilterMemory[column][COHPSK_NSYM - 1] = ComplexMath.times(baseband[column], GAIN);    // bottom row
        }

        /*
         * Run the audio filter over the constellation data for each row
         */
        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            for (int i = 0; i < COHPSK_M; i++) {                // 0..99
                Complex acc = new Complex();

                for (int j = 0, k = COHPSK_M - i - 1; j < COHPSK_NSYM; j++, k += COHPSK_M) {
                    acc = ComplexMath.add(acc, ComplexMath.times(ComplexMath.times(txFilterMemory[column][j], COHPSK_M), GTALPHA5ROOT[k]));
                }

                // Adjust the baseband phase and add symbol
                mtxPhase[column] = ComplexMath.times(mtxPhase[column], mtxFrequency[column]);
                waveform[offset + i] = ComplexMath.add(waveform[offset + i], ComplexMath.times(acc, mtxPhase[column]));
            }
        }

        // Adjust the final phase and move the baseband segment up to 1500 Hz
        for (int i = 0; i < COHPSK_M; i++) {
            mbasebandTxPhase = ComplexMath.times(mbasebandTxPhase, TXPHASE);
            waveform[offset + i] = ComplexMath.times(ComplexMath.times(waveform[offset + i], mbasebandTxPhase), TWO);
        }

        // Normalize carriers of the baseband phase for drift
        for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
            mtxPhase[column] = ComplexMath.normalize(mtxPhase[column]);
        }

        // same for upconversion phase, combat drift
        mbasebandTxPhase = ComplexMath.normalize(mbasebandTxPhase);

        /*
         * Move the rows 0..4 up one older position (in time). You don't
         * have to worry about zero'ing out the last row, as the new rows
         * carrier symbols will go there on the next call to this method.
         */
        for (int row = 0; row < (COHPSK_NSYM - 1); row++) {
            for (int column = 0; column < (COHPSK_NC * mdiversityFactor); column++) {
                txFilterMemory[column][row] = txFilterMemory[column][row + 1];
            }
        }
    }
}
