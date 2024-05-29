/*
 * Copyright (C) 2015 Jim Ahlstrom
 *
 * All rights reserved.
 *
 * Java Version by Steve Sampson
 */
package cohfdm;

import complex.Complex;
import complex.ComplexMath;

/**
 * Sample rate converter
 *
 * This filter is from the Codec2 project by David Rowe, but donated to the
 * project by Jim Ahlstrom from his Quisk SDR project.
 *
 * @author http://james.ahlstrom.name/quisk/
 */
public final class Resample implements IDefines {

    private final Complex[] mSamples;
    private final Complex[] mBuffer;
    private int msamplePointer;
    private int mdecimationIndex;
    private final int mfilterTaps;

    public Resample() {
        mdecimationIndex = 0;   // index of next sample for decimation
        msamplePointer = 0;     // next available position in mSamples buffer

        mBuffer = new Complex[COHPSK_NOM_RX_SAMPLES_PER_FRAME];  // buffer for interpolation

        for (int i = 0; i < mBuffer.length; i++) {
            mBuffer[i] = new Complex();
        }

        mfilterTaps = COEFS.length;
        mSamples = new Complex[mfilterTaps];

        for (int i = 0; i < mSamples.length; i++) {
            mSamples[i] = new Complex();
        }
    }

    /**
     * Take an array of complex Samples of length count, multiply the sample
     * rate by interp, and then divide the sample rate by decim.
     *
     * Return the new number of samples. Each specific interp and decim will
     * require its own custom FIR filter.
     *
     * @param cSamples complex modem samples
     * @param count int of the number of samples to use
     * @param interp int interpolation factor
     * @param decim int decimation factor
     * @return int of the resulting number of samples
     */
    public int interpDecim(Complex[] cSamples, int count, int interp, int decim) {
        // Copy cSamples to auxillary buffer, because we are
        // going to re-use cSamples for output
        System.arraycopy(cSamples, 0, mBuffer, 0, count);

        int nOut = 0;

        for (int i = 0; i < count; i++) {
            mSamples[msamplePointer] = mBuffer[i];

            while (mdecimationIndex < interp) {
                int ptSample = msamplePointer;
                int ptCoef = mdecimationIndex;
                Complex csample = new Complex();// complex class is immutable
                // so can't just reset to zero.
                for (int k = 0; k < (mfilterTaps / interp); k++, ptCoef += interp) {
                    csample = ComplexMath.add(csample,
                            ComplexMath.times(mSamples[ptSample], COEFS[ptCoef])
                    );

                    ptSample--;

                    if (ptSample < 0) {
                        ptSample = mfilterTaps - 1;
                    }
                }

                cSamples[nOut] = ComplexMath.times(csample, interp);
                nOut++;
                mdecimationIndex += decim;
            }

            msamplePointer++;

            if (msamplePointer >= mfilterTaps) {
                msamplePointer = 0;
            }

            mdecimationIndex -= interp;
        }

        return nOut;
    }

    /**
     * Getter for Index of next sample for decimation
     *
     * @return int next decimation index
     */
    public int getDecimationIndex() {
        return mdecimationIndex;
    }
}
