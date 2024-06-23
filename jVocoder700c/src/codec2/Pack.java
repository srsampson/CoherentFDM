/*
 * Copyright (C) 2010 Perens LLC
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 *
 * Java Version by Steve Sampson
 */
package codec2;

public final class Pack {

    private int m_bitOffset;

    protected Pack() {
        m_bitOffset = 0;
    }

    protected void reset() {
        m_bitOffset = 0;
    }

    protected void pack(byte[] bitArray, int value, int bits) {
        int valueBits = bits;

        do {
            int bI = m_bitOffset;
            int bitsLeft = 8 - (bI & 0x7);
            int sliceWidth = bitsLeft < valueBits ? bitsLeft : valueBits;
            int wordIndex = bI >>> 3;

            bitArray[wordIndex] |= ((byte) ((value >> (valueBits - sliceWidth)) << (bitsLeft - sliceWidth)));

            m_bitOffset = bI + sliceWidth;
            valueBits -= sliceWidth;
        } while (valueBits != 0);
    }

    protected int unpack(byte[] bitArray, int bits) {
        int valueBits = bits;
        int field = 0;

        do {
            int bI = m_bitOffset;
            int bitsLeft = 8 - (bI & 0x7);
            int sliceWidth = bitsLeft < valueBits ? bitsLeft : valueBits;

            field |= (((bitArray[bI >>> 3] >> (bitsLeft - sliceWidth)) & ((1 << sliceWidth) - 1)) << (valueBits - sliceWidth));

            m_bitOffset = bI + sliceWidth;
            valueBits -= sliceWidth;
        } while (valueBits != 0);

        return field;
    }
}
