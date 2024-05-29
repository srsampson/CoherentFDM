/*
 * Copyright (C) 2007-2009 Piotr Wendykier, Emory University
 *
 * All Rights Reserved
 *
 * The Original Code is JTransforms.
 *
 * This code is derived from General Purpose FFT Package written by Takuya Ooura
 * and from JFFTPack written by Baoshe Zhang
 */
package codec2;

import complex.Complex;

public final class FFT implements IDefines {

    public static final int[] FACTORS = {4, 2, 3, 5};
    //
    private final float[] m_wtable;
    private final float[] m_wtable_r;
    private final float[] m_c;
    private final int m_n;

    public FFT(int size) {
        m_n = size;
        m_wtable = new float[4 * m_n + 15];
        m_wtable_r = new float[2 * m_n + 15];
        m_c = new float[m_n * 2];

        cffti();
        rffti();
    }

    public void transform(Complex[] a) {
        for (int i = 0, j = 0; i < (m_n * 2); i += 2, j++) {
            m_c[i] = a[j].real();
            m_c[i + 1] = a[j].imag();
        }

        cfftf(m_c, -1);

        for (int i = 0, j = 0; i < (m_n * 2); i += 2, j++) {
            a[j] = new Complex(m_c[i], m_c[i + 1]);
        }
    }

    public void itransform(Complex[] a) {
        for (int i = 0, j = 0; i < (m_n * 2); i += 2, j++) {
            m_c[i] = a[j].real();
            m_c[i + 1] = a[j].imag();
        }

        cfftf(m_c, +1);

        for (int i = 0, j = 0; i < (m_n * 2); i += 2, j++) {
            a[j] = new Complex(m_c[i], m_c[i + 1]);
        }
    }

   /*---------------------------------------------------------
     cffti: initialization of Complex FFT
     --------------------------------------------------------*/
    private void cffti() {
        boolean loop;
        int twon = 2 * m_n;
        int fourn = 4 * m_n;

        int nl = m_n;
        int nf = 0;
        int j = 0;
        int ntry = 0;

        do {
            loop = false;

            j++;

            if (j <= 4) {
                ntry = FACTORS[j - 1];
            } else {
                ntry += 2;
            }

            do {
                int nq = nl / ntry;
                int nr = nl - ntry * nq;

                if (nr != 0) {
                    loop = true;
                    break;
                }

                nf++;
                m_wtable[nf + 1 + fourn] = ntry;
                nl = nq;

                if (ntry == 2 && nf != 1) {
                    for (int i = 2; i <= nf; i++) {
                        int ib = nf - i + 2;
                        int idx = ib + fourn;
                        m_wtable[idx + 1] = m_wtable[idx];
                    }

                    m_wtable[2 + fourn] = 2;
                }
            } while (nl != 1);
        } while (loop);

        m_wtable[fourn] = m_n;
        m_wtable[1 + fourn] = nf;

        float argh = TAU / (float) m_n;
        int iz = 1;
        int l1 = 1;

        for (int k1 = 1; k1 <= nf; k1++) {
            int ip = (int) m_wtable[k1 + 1 + fourn];
            int ld = 0;
            int l2 = l1 * ip;

            int ido = m_n / l2;
            int idot = ido + ido + 2;
            int ipm = ip - 1;

            for (j = 1; j <= ipm; j++) {
                int i1 = iz;
                m_wtable[iz - 1 + twon] = 1;
                m_wtable[iz + twon] = 0;

                ld += l1;
                float fi = 0.0f;
                float argld = ld * argh;

                for (int ii = 4; ii <= idot; ii += 2) {
                    iz += 2;
                    fi += 1.0f;
                    float arg = fi * argld;
                    int idx = iz + twon;
                    m_wtable[idx - 1] = (float) Math.cos(arg);
                    m_wtable[idx] = (float) Math.sin(arg);
                }

                if (ip > 5) {
                    int idx1 = i1 + twon;
                    int idx2 = iz + twon;
                    m_wtable[idx1 - 1] = m_wtable[idx2 - 1];
                    m_wtable[idx1] = m_wtable[idx2];
                }
            }

            l1 = l2;
        }
    }

    private void rffti() {
        boolean loop;
        int twon = 2 * m_n;
        int nl = m_n;
        int nf = 0;
        int j = 0;
        int ntry = 0;

        do {
            loop = false;

            ++j;

            if (j <= 4) {
                ntry = FACTORS[j - 1];
            } else {
                ntry += 2;
            }

            do {
                int nq = nl / ntry;
                int nr = nl - ntry * nq;

                if (nr != 0) {
                    loop = true;
                    break;
                }

                ++nf;
                m_wtable_r[nf + 1 + twon] = ntry;

                nl = nq;
                if (ntry == 2 && nf != 1) {
                    for (int i = 2; i <= nf; i++) {
                        int ib = nf - i + 2;
                        int idx = ib + twon;
                        m_wtable_r[idx + 1] = m_wtable_r[idx];
                    }

                    m_wtable_r[2 + twon] = 2;
                }
            } while (nl != 1);
        } while (loop);

        m_wtable_r[twon] = m_n;
        m_wtable_r[1 + twon] = nf;

        float argh = TAU / (float) m_n;
        int is = 0;

        int nfm1 = nf - 1;
        int l1 = 1;

        if (nfm1 == 0) {
            return;
        }

        for (int k1 = 1; k1 <= nfm1; k1++) {
            int ip = (int) m_wtable_r[k1 + 1 + twon];
            int ld = 0;
            int l2 = l1 * ip;
            int ido = m_n / l2;
            int ipm = ip - 1;

            for (j = 1; j <= ipm; ++j) {
                ld += l1;
                int iz = is;
                float argld = ld * argh;
                float fi = 0.0f;

                for (int ii = 3; ii <= ido; ii += 2) {
                    iz += 2;
                    fi += 1;
                    float arg = fi * argld;
                    int idx = iz + m_n;
                    m_wtable_r[idx - 2] = (float) Math.cos(arg);
                    m_wtable_r[idx - 1] = (float) Math.sin(arg);
                }

                is += ido;
            }

            l1 = l2;
        }
    }

    /*---------------------------------------------------------
     cfftf1: further processing of Complex forward FFT
     --------------------------------------------------------*/
    private void cfftf(float[] a, int isign) {
        boolean[] nac = new boolean[1];

        int twon = 2 * m_n;
        float[] ch = new float[twon];

        nac[0] = false;

        int iw1 = twon;
        int iw = iw1;

        int iw2 = 4 * m_n;
        int nf = (int) m_wtable[1 + iw2];

        int na = 0;
        int l1 = 1;

        for (int k1 = 2; k1 <= nf + 1; k1++) {
            int ip = (int) m_wtable[k1 + iw2];
            int l2 = ip * l1;
            int ido = m_n / l2;
            int idot = ido + ido;
            int idl1 = idot * l1;

            switch (ip) {
                case 4:
                    if (na == 0) {
                        passf4(idot, l1, a, 0, ch, 0, iw, isign);
                    } else {
                        passf4(idot, l1, ch, 0, a, 0, iw, isign);
                    }
                    na = 1 - na;
                    break;
                case 2:
                    if (na == 0) {
                        passf2(idot, l1, a, 0, ch, 0, iw, isign);
                    } else {
                        passf2(idot, l1, ch, 0, a, 0, iw, isign);
                    }
                    na = 1 - na;
                    break;
                case 3:
                    if (na == 0) {
                        passf3(idot, l1, a, 0, ch, 0, iw, isign);
                    } else {
                        passf3(idot, l1, ch, 0, a, 0, iw, isign);
                    }
                    na = 1 - na;
                    break;
                case 5:
                    if (na == 0) {
                        passf5(idot, l1, a, 0, ch, 0, iw, isign);
                    } else {
                        passf5(idot, l1, ch, 0, a, 0, iw, isign);
                    }
                    na = 1 - na;
                    break;
                default:
                    if (na == 0) {
                        passfg(nac, idot, ip, l1, idl1, a, 0, ch, 0, iw, isign);
                    } else {
                        passfg(nac, idot, ip, l1, idl1, ch, 0, a, 0, iw, isign);
                    }
                    if (nac[0] == true) {
                        na = 1 - na;
                    }
                    break;
            }

            l1 = l2;
            iw += (ip - 1) * idot;
        }

        if (na == 0) {
            return;
        }

        System.arraycopy(ch, 0, a, 0, twon);
    }

    /*----------------------------------------------------------------------
     passf2: Complex FFT's forward/backward processing of factor 2;
     isign is +1 for backward and -1 for forward transforms
     ----------------------------------------------------------------------*/
    private void passf2(int ido, int l1, float[] in, int in_off, float[] out, int out_off, int offset, int isign) {
        int iw1 = offset;
        int idx = ido * l1;

        if (ido <= 2) {
            for (int k = 0; k < l1; k++) {
                int idx0 = k * ido;
                int iidx1 = in_off + 2 * idx0;
                int iidx2 = iidx1 + ido;
                float a1r = in[iidx1];
                float a1i = in[iidx1 + 1];
                float a2r = in[iidx2];
                float a2i = in[iidx2 + 1];

                int oidx1 = out_off + idx0;
                int oidx2 = oidx1 + idx;
                out[oidx1] = a1r + a2r;
                out[oidx1 + 1] = a1i + a2i;
                out[oidx2] = a1r - a2r;
                out[oidx2 + 1] = a1i - a2i;
            }
        } else {
            for (int k = 0; k < l1; k++) {
                for (int i = 0; i < ido - 1; i += 2) {
                    int idx0 = k * ido;
                    int iidx1 = in_off + i + 2 * idx0;
                    int iidx2 = iidx1 + ido;
                    float i1r = in[iidx1];
                    float i1i = in[iidx1 + 1];
                    float i2r = in[iidx2];
                    float i2i = in[iidx2 + 1];

                    int widx1 = i + iw1;
                    float w1r = m_wtable[widx1];
                    float w1i = isign * m_wtable[widx1 + 1];

                    float t1r = i1r - i2r;
                    float t1i = i1i - i2i;

                    int oidx1 = out_off + i + idx0;
                    int oidx2 = oidx1 + idx;
                    out[oidx1] = i1r + i2r;
                    out[oidx1 + 1] = i1i + i2i;
                    out[oidx2] = w1r * t1r - w1i * t1i;
                    out[oidx2 + 1] = w1r * t1i + w1i * t1r;
                }
            }
        }
    }

    /*----------------------------------------------------------------------
     passf3: Complex FFT's forward/backward processing of factor 3;
     isign is +1 for backward and -1 for forward transforms
     ----------------------------------------------------------------------*/
    private void passf3(int ido, int l1, float[] in, int in_off, float[] out, int out_off, int offset, int isign) {
        float taur = -0.5f;
        float taui = 0.866025403784438707610604524234076962f;
        float ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

        int iw1 = offset;
        int iw2 = iw1 + ido;
        int idxt = l1 * ido;

        if (ido == 2) {
            for (int k = 1; k <= l1; k++) {
                int iidx1 = in_off + (3 * k - 2) * ido;
                int iidx2 = iidx1 + ido;
                int iidx3 = iidx1 - ido;
                float i1r = in[iidx1];
                float i1i = in[iidx1 + 1];
                float i2r = in[iidx2];
                float i2i = in[iidx2 + 1];
                float i3r = in[iidx3];
                float i3i = in[iidx3 + 1];

                tr2 = i1r + i2r;
                cr2 = i3r + taur * tr2;
                ti2 = i1i + i2i;
                ci2 = i3i + taur * ti2;
                cr3 = isign * taui * (i1r - i2r);
                ci3 = isign * taui * (i1i - i2i);

                int oidx1 = out_off + (k - 1) * ido;
                int oidx2 = oidx1 + idxt;
                int oidx3 = oidx2 + idxt;
                out[oidx1] = in[iidx3] + tr2;
                out[oidx1 + 1] = i3i + ti2;
                out[oidx2] = cr2 - ci3;
                out[oidx2 + 1] = ci2 + cr3;
                out[oidx3] = cr2 + ci3;
                out[oidx3 + 1] = ci2 - cr3;
            }
        } else {
            for (int k = 1; k <= l1; k++) {
                int idx1 = in_off + (3 * k - 2) * ido;
                int idx2 = out_off + (k - 1) * ido;
                for (int i = 0; i < ido - 1; i += 2) {
                    int iidx1 = i + idx1;
                    int iidx2 = iidx1 + ido;
                    int iidx3 = iidx1 - ido;
                    float a1r = in[iidx1];
                    float a1i = in[iidx1 + 1];
                    float a2r = in[iidx2];
                    float a2i = in[iidx2 + 1];
                    float a3r = in[iidx3];
                    float a3i = in[iidx3 + 1];

                    tr2 = a1r + a2r;
                    cr2 = a3r + taur * tr2;
                    ti2 = a1i + a2i;
                    ci2 = a3i + taur * ti2;
                    cr3 = isign * taui * (a1r - a2r);
                    ci3 = isign * taui * (a1i - a2i);
                    dr2 = cr2 - ci3;
                    dr3 = cr2 + ci3;
                    di2 = ci2 + cr3;
                    di3 = ci2 - cr3;

                    int widx1 = i + iw1;
                    int widx2 = i + iw2;
                    float w1r = m_wtable[widx1];
                    float w1i = isign * m_wtable[widx1 + 1];
                    float w2r = m_wtable[widx2];
                    float w2i = isign * m_wtable[widx2 + 1];

                    int oidx1 = i + idx2;
                    int oidx2 = oidx1 + idxt;
                    int oidx3 = oidx2 + idxt;
                    out[oidx1] = a3r + tr2;
                    out[oidx1 + 1] = a3i + ti2;
                    out[oidx2] = w1r * dr2 - w1i * di2;
                    out[oidx2 + 1] = w1r * di2 + w1i * dr2;
                    out[oidx3] = w2r * dr3 - w2i * di3;
                    out[oidx3 + 1] = w2r * di3 + w2i * dr3;
                }
            }
        }
    }

    /*----------------------------------------------------------------------
     passf4: Complex FFT's forward/backward processing of factor 4;
     isign is +1 for backward and -1 for forward transforms
     ----------------------------------------------------------------------*/
    private void passf4(int ido, int l1, float[] in, int in_off, float[] out, int out_off, int offset, int isign) {
        float ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
        int iw1 = offset;
        int iw2 = iw1 + ido;
        int iw3 = iw2 + ido;
        int idx0 = l1 * ido;

        if (ido == 2) {
            for (int k = 0; k < l1; k++) {
                int idxt1 = k * ido;
                int iidx1 = in_off + 4 * idxt1 + 1;
                int iidx2 = iidx1 + ido;
                int iidx3 = iidx2 + ido;
                int iidx4 = iidx3 + ido;

                float i1i = in[iidx1 - 1];
                float i1r = in[iidx1];
                float i2i = in[iidx2 - 1];
                float i2r = in[iidx2];
                float i3i = in[iidx3 - 1];
                float i3r = in[iidx3];
                float i4i = in[iidx4 - 1];
                float i4r = in[iidx4];

                ti1 = i1r - i3r;
                ti2 = i1r + i3r;
                tr4 = i4r - i2r;
                ti3 = i2r + i4r;
                tr1 = i1i - i3i;
                tr2 = i1i + i3i;
                ti4 = i2i - i4i;
                tr3 = i2i + i4i;

                int oidx1 = out_off + idxt1;
                int oidx2 = oidx1 + idx0;
                int oidx3 = oidx2 + idx0;
                int oidx4 = oidx3 + idx0;
                out[oidx1] = tr2 + tr3;
                out[oidx1 + 1] = ti2 + ti3;
                out[oidx2] = tr1 + isign * tr4;
                out[oidx2 + 1] = ti1 + isign * ti4;
                out[oidx3] = tr2 - tr3;
                out[oidx3 + 1] = ti2 - ti3;
                out[oidx4] = tr1 - isign * tr4;
                out[oidx4 + 1] = ti1 - isign * ti4;
            }
        } else {
            for (int k = 0; k < l1; k++) {
                int idx1 = k * ido;
                int idx2 = in_off + 1 + 4 * idx1;
                for (int i = 0; i < ido - 1; i += 2) {
                    int iidx1 = i + idx2;
                    int iidx2 = iidx1 + ido;
                    int iidx3 = iidx2 + ido;
                    int iidx4 = iidx3 + ido;
                    float i1i = in[iidx1 - 1];
                    float i1r = in[iidx1];
                    float i2i = in[iidx2 - 1];
                    float i2r = in[iidx2];
                    float i3i = in[iidx3 - 1];
                    float i3r = in[iidx3];
                    float i4i = in[iidx4 - 1];
                    float i4r = in[iidx4];

                    ti1 = i1r - i3r;
                    ti2 = i1r + i3r;
                    ti3 = i2r + i4r;
                    tr4 = i4r - i2r;
                    tr1 = i1i - i3i;
                    tr2 = i1i + i3i;
                    ti4 = i2i - i4i;
                    tr3 = i2i + i4i;
                    cr3 = tr2 - tr3;
                    ci3 = ti2 - ti3;
                    cr2 = tr1 + isign * tr4;
                    cr4 = tr1 - isign * tr4;
                    ci2 = ti1 + isign * ti4;
                    ci4 = ti1 - isign * ti4;

                    int widx1 = i + iw1;
                    int widx2 = i + iw2;
                    int widx3 = i + iw3;
                    float w1r = m_wtable[widx1];
                    float w1i = isign * m_wtable[widx1 + 1];
                    float w2r = m_wtable[widx2];
                    float w2i = isign * m_wtable[widx2 + 1];
                    float w3r = m_wtable[widx3];
                    float w3i = isign * m_wtable[widx3 + 1];

                    int oidx1 = out_off + i + idx1;
                    int oidx2 = oidx1 + idx0;
                    int oidx3 = oidx2 + idx0;
                    int oidx4 = oidx3 + idx0;
                    out[oidx1] = tr2 + tr3;
                    out[oidx1 + 1] = ti2 + ti3;
                    out[oidx2] = w1r * cr2 - w1i * ci2;
                    out[oidx2 + 1] = w1r * ci2 + w1i * cr2;
                    out[oidx3] = w2r * cr3 - w2i * ci3;
                    out[oidx3 + 1] = w2r * ci3 + w2i * cr3;
                    out[oidx4] = w3r * cr4 - w3i * ci4;
                    out[oidx4 + 1] = w3r * ci4 + w3i * cr4;
                }
            }
        }
    }

    /*----------------------------------------------------------------------
     passf5: Complex FFT's forward/backward processing of factor 5;
     isign is +1 for backward and -1 for forward transforms
     ----------------------------------------------------------------------*/
    private void passf5(int ido, int l1, float[] in, int in_off, float[] out, int out_off, int offset, int isign) /* isign==-1 for forward transform and+1 for backward transform */ {
        float tr11 = 0.309016994374947451262869435595348477f;
        float ti11 = 0.951056516295153531181938433292089030f;
        float tr12 = -0.809016994374947340240566973079694435f;
        float ti12 = 0.587785252292473248125759255344746634f;
        float ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

        int iw1 = offset;
        int iw2 = iw1 + ido;
        int iw3 = iw2 + ido;
        int iw4 = iw3 + ido;
        int idx0 = l1 * ido;

        if (ido == 2) {
            for (int k = 1; k <= l1; ++k) {
                int iidx1 = in_off + (5 * k - 4) * ido + 1;
                int iidx2 = iidx1 + ido;
                int iidx3 = iidx1 - ido;
                int iidx4 = iidx2 + ido;
                int iidx5 = iidx4 + ido;

                float i1i = in[iidx1 - 1];
                float i1r = in[iidx1];
                float i2i = in[iidx2 - 1];
                float i2r = in[iidx2];
                float i3i = in[iidx3 - 1];
                float i3r = in[iidx3];
                float i4i = in[iidx4 - 1];
                float i4r = in[iidx4];
                float i5i = in[iidx5 - 1];
                float i5r = in[iidx5];

                ti5 = i1r - i5r;
                ti2 = i1r + i5r;
                ti4 = i2r - i4r;
                ti3 = i2r + i4r;
                tr5 = i1i - i5i;
                tr2 = i1i + i5i;
                tr4 = i2i - i4i;
                tr3 = i2i + i4i;
                cr2 = i3i + tr11 * tr2 + tr12 * tr3;
                ci2 = i3r + tr11 * ti2 + tr12 * ti3;
                cr3 = i3i + tr12 * tr2 + tr11 * tr3;
                ci3 = i3r + tr12 * ti2 + tr11 * ti3;
                cr5 = isign * (ti11 * tr5 + ti12 * tr4);
                ci5 = isign * (ti11 * ti5 + ti12 * ti4);
                cr4 = isign * (ti12 * tr5 - ti11 * tr4);
                ci4 = isign * (ti12 * ti5 - ti11 * ti4);

                int oidx1 = out_off + (k - 1) * ido;
                int oidx2 = oidx1 + idx0;
                int oidx3 = oidx2 + idx0;
                int oidx4 = oidx3 + idx0;
                int oidx5 = oidx4 + idx0;
                out[oidx1] = i3i + tr2 + tr3;
                out[oidx1 + 1] = i3r + ti2 + ti3;
                out[oidx2] = cr2 - ci5;
                out[oidx2 + 1] = ci2 + cr5;
                out[oidx3] = cr3 - ci4;
                out[oidx3 + 1] = ci3 + cr4;
                out[oidx4] = cr3 + ci4;
                out[oidx4 + 1] = ci3 - cr4;
                out[oidx5] = cr2 + ci5;
                out[oidx5 + 1] = ci2 - cr5;
            }
        } else {
            for (int k = 1; k <= l1; k++) {
                int idx1 = in_off + 1 + (k * 5 - 4) * ido;
                int idx2 = out_off + (k - 1) * ido;
                for (int i = 0; i < ido - 1; i += 2) {
                    int iidx1 = i + idx1;
                    int iidx2 = iidx1 + ido;
                    int iidx3 = iidx1 - ido;
                    int iidx4 = iidx2 + ido;
                    int iidx5 = iidx4 + ido;
                    float i1i = in[iidx1 - 1];
                    float i1r = in[iidx1];
                    float i2i = in[iidx2 - 1];
                    float i2r = in[iidx2];
                    float i3i = in[iidx3 - 1];
                    float i3r = in[iidx3];
                    float i4i = in[iidx4 - 1];
                    float i4r = in[iidx4];
                    float i5i = in[iidx5 - 1];
                    float i5r = in[iidx5];

                    ti5 = i1r - i5r;
                    ti2 = i1r + i5r;
                    ti4 = i2r - i4r;
                    ti3 = i2r + i4r;
                    tr5 = i1i - i5i;
                    tr2 = i1i + i5i;
                    tr4 = i2i - i4i;
                    tr3 = i2i + i4i;
                    cr2 = i3i + tr11 * tr2 + tr12 * tr3;
                    ci2 = i3r + tr11 * ti2 + tr12 * ti3;
                    cr3 = i3i + tr12 * tr2 + tr11 * tr3;
                    ci3 = i3r + tr12 * ti2 + tr11 * ti3;
                    cr5 = isign * (ti11 * tr5 + ti12 * tr4);
                    ci5 = isign * (ti11 * ti5 + ti12 * ti4);
                    cr4 = isign * (ti12 * tr5 - ti11 * tr4);
                    ci4 = isign * (ti12 * ti5 - ti11 * ti4);
                    dr3 = cr3 - ci4;
                    dr4 = cr3 + ci4;
                    di3 = ci3 + cr4;
                    di4 = ci3 - cr4;
                    dr5 = cr2 + ci5;
                    dr2 = cr2 - ci5;
                    di5 = ci2 - cr5;
                    di2 = ci2 + cr5;

                    int widx1 = i + iw1;
                    int widx2 = i + iw2;
                    int widx3 = i + iw3;
                    int widx4 = i + iw4;
                    float w1r = m_wtable[widx1];
                    float w1i = isign * m_wtable[widx1 + 1];
                    float w2r = m_wtable[widx2];
                    float w2i = isign * m_wtable[widx2 + 1];
                    float w3r = m_wtable[widx3];
                    float w3i = isign * m_wtable[widx3 + 1];
                    float w4r = m_wtable[widx4];
                    float w4i = isign * m_wtable[widx4 + 1];

                    int oidx1 = i + idx2;
                    int oidx2 = oidx1 + idx0;
                    int oidx3 = oidx2 + idx0;
                    int oidx4 = oidx3 + idx0;
                    int oidx5 = oidx4 + idx0;
                    out[oidx1] = i3i + tr2 + tr3;
                    out[oidx1 + 1] = i3r + ti2 + ti3;
                    out[oidx2] = w1r * dr2 - w1i * di2;
                    out[oidx2 + 1] = w1r * di2 + w1i * dr2;
                    out[oidx3] = w2r * dr3 - w2i * di3;
                    out[oidx3 + 1] = w2r * di3 + w2i * dr3;
                    out[oidx4] = w3r * dr4 - w3i * di4;
                    out[oidx4 + 1] = w3r * di4 + w3i * dr4;
                    out[oidx5] = w4r * dr5 - w4i * di5;
                    out[oidx5 + 1] = w4r * di5 + w4i * dr5;
                }
            }
        }
    }

    /*----------------------------------------------------------------------
     passfg: Complex FFT's forward/backward processing of general factor;
     isign is +1 for backward and -1 for forward transforms
     ----------------------------------------------------------------------*/
    private void passfg(boolean[] nac, int ido, int ip, int l1, int idl1, float[] in, int in_off, float[] out, int out_off, int offset, int isign) {
        int idij, idlj, l, jc, lc, idj;
        float w1r, w1i, w2i, w2r;

        int iw1 = offset;
        int idot = ido / 2;
        int ipph = (ip + 1) / 2;
        int idp = ip * ido;

        if (ido >= l1) {
            for (int j = 1; j < ipph; j++) {
                jc = ip - j;
                int idx1 = j * ido;
                int idx2 = jc * ido;

                for (int k = 0; k < l1; k++) {
                    int idx3 = k * ido;
                    int idx4 = idx3 + idx1 * l1;
                    int idx5 = idx3 + idx2 * l1;
                    int idx6 = idx3 * ip;

                    for (int i = 0; i < ido; i++) {
                        int oidx1 = out_off + i;
                        float i1r = in[in_off + i + idx1 + idx6];
                        float i2r = in[in_off + i + idx2 + idx6];
                        out[oidx1 + idx4] = i1r + i2r;
                        out[oidx1 + idx5] = i1r - i2r;
                    }
                }
            }

            for (int k = 0; k < l1; k++) {
                int idxt1 = k * ido;
                int idxt2 = idxt1 * ip;

                for (int i = 0; i < ido; i++) {
                    out[out_off + i + idxt1] = in[in_off + i + idxt2];
                }
            }
        } else {
            for (int j = 1; j < ipph; j++) {
                jc = ip - j;
                int idxt1 = j * l1 * ido;
                int idxt2 = jc * l1 * ido;
                int idxt3 = j * ido;
                int idxt4 = jc * ido;

                for (int i = 0; i < ido; i++) {
                    for (int k = 0; k < l1; k++) {
                        int idx1 = k * ido;
                        int idx2 = idx1 * ip;
                        int idx3 = out_off + i;
                        int idx4 = in_off + i;
                        float i1r = in[idx4 + idxt3 + idx2];
                        float i2r = in[idx4 + idxt4 + idx2];
                        out[idx3 + idx1 + idxt1] = i1r + i2r;
                        out[idx3 + idx1 + idxt2] = i1r - i2r;
                    }
                }
            }

            for (int i = 0; i < ido; i++) {
                for (int k = 0; k < l1; k++) {
                    int idx1 = k * ido;
                    out[out_off + i + idx1] = in[in_off + i + idx1 * ip];
                }
            }
        }

        int idl = 2 - ido;
        int inc = 0;
        int idxt0 = (ip - 1) * idl1;

        for (l = 1; l < ipph; l++) {
            lc = ip - l;
            idl += ido;
            int idxt1 = l * idl1;
            int idxt2 = lc * idl1;
            int idxt3 = idl + iw1;
            w1r = m_wtable[idxt3 - 2];
            w1i = isign * m_wtable[idxt3 - 1];

            for (int ik = 0; ik < idl1; ik++) {
                int idx1 = in_off + ik;
                int idx2 = out_off + ik;
                in[idx1 + idxt1] = out[idx2] + w1r * out[idx2 + idl1];
                in[idx1 + idxt2] = w1i * out[idx2 + idxt0];
            }

            idlj = idl;
            inc += ido;

            for (int j = 2; j < ipph; j++) {
                jc = ip - j;
                idlj += inc;

                if (idlj > idp) {
                    idlj -= idp;
                }

                int idxt4 = idlj + iw1;
                w2r = m_wtable[idxt4 - 2];
                w2i = isign * m_wtable[idxt4 - 1];
                int idxt5 = j * idl1;
                int idxt6 = jc * idl1;

                for (int ik = 0; ik < idl1; ik++) {
                    int idx1 = in_off + ik;
                    int idx2 = out_off + ik;
                    in[idx1 + idxt1] += w2r * out[idx2 + idxt5];
                    in[idx1 + idxt2] += w2i * out[idx2 + idxt6];
                }
            }
        }

        for (int j = 1; j < ipph; j++) {
            int idxt1 = j * idl1;
            for (int ik = 0; ik < idl1; ik++) {
                int idx1 = out_off + ik;
                out[idx1] += out[idx1 + idxt1];
            }
        }

        for (int j = 1; j < ipph; j++) {
            jc = ip - j;
            int idx1 = j * idl1;
            int idx2 = jc * idl1;

            for (int ik = 1; ik < idl1; ik += 2) {
                int idx3 = out_off + ik;
                int idx4 = in_off + ik;
                int iidx1 = idx4 + idx1;
                int iidx2 = idx4 + idx2;
                float i1i = in[iidx1 - 1];
                float i1r = in[iidx1];
                float i2i = in[iidx2 - 1];
                float i2r = in[iidx2];

                int oidx1 = idx3 + idx1;
                int oidx2 = idx3 + idx2;
                out[oidx1 - 1] = i1i - i2r;
                out[oidx2 - 1] = i1i + i2r;
                out[oidx1] = i1r + i2i;
                out[oidx2] = i1r - i2i;
            }
        }

        nac[0] = true;

        if (ido == 2) {
            return;
        }

        nac[0] = false;
        System.arraycopy(out, out_off, in, in_off, idl1);
        int idx0 = l1 * ido;

        for (int j = 1; j < ip; j++) {
            int idx1 = j * idx0;

            for (int k = 0; k < l1; k++) {
                int idx2 = k * ido;
                int oidx1 = out_off + idx2 + idx1;
                int iidx1 = in_off + idx2 + idx1;
                in[iidx1] = out[oidx1];
                in[iidx1 + 1] = out[oidx1 + 1];
            }
        }

        if (idot <= l1) {
            idij = 0;

            for (int j = 1; j < ip; j++) {
                idij += 2;
                int idx1 = j * l1 * ido;

                for (int i = 3; i < ido; i += 2) {
                    idij += 2;
                    int idx2 = idij + iw1 - 1;
                    w1r = m_wtable[idx2 - 1];
                    w1i = isign * m_wtable[idx2];
                    int idx3 = in_off + i;
                    int idx4 = out_off + i;

                    for (int k = 0; k < l1; k++) {
                        int idx5 = k * ido + idx1;
                        int iidx1 = idx3 + idx5;
                        int oidx1 = idx4 + idx5;
                        float o1i = out[oidx1 - 1];
                        float o1r = out[oidx1];
                        in[iidx1 - 1] = w1r * o1i - w1i * o1r;
                        in[iidx1] = w1r * o1r + w1i * o1i;
                    }
                }
            }
        } else {
            idj = 2 - ido;

            for (int j = 1; j < ip; j++) {
                idj += ido;
                int idx1 = j * l1 * ido;

                for (int k = 0; k < l1; k++) {
                    idij = idj;
                    int idx3 = k * ido + idx1;

                    for (int i = 3; i < ido; i += 2) {
                        idij += 2;
                        int idx2 = idij - 1 + iw1;
                        w1r = m_wtable[idx2 - 1];
                        w1i = isign * m_wtable[idx2];
                        int iidx1 = in_off + i + idx3;
                        int oidx1 = out_off + i + idx3;
                        float o1i = out[oidx1 - 1];
                        float o1r = out[oidx1];
                        in[iidx1 - 1] = w1r * o1i - w1i * o1r;
                        in[iidx1] = w1r * o1r + w1i * o1i;
                    }
                }
            }
        }
    }
}
