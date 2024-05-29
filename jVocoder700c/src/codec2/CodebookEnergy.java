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

public final class CodebookEnergy {

    private static final float[] CODES0 = {
        10.0f, 12.5f, 15.0f, 17.5f, 20.0f,
        22.5f, 25.0f, 27.5f, 30.0f, 32.5f,
        35.0f, 37.5f, 40.0f, 42.5f, 45.0f,
        47.5f
    };

    private final Codebook m_energy_cb;

    public CodebookEnergy() {
        m_energy_cb = new Codebook(1, 4, 16, CODES0);
    }

    public Codebook getCodebook() {
        return m_energy_cb;
    }
}
