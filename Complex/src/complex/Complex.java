/*
 * Java by Steve Sampson, K5OKC
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */
package complex;

/**
 * @author Steve Sampson, K5OKC
 * @version .1-pre-alpha, 03/18/2023
 */
public class Complex {

    private final float real;
    private final float imag;

    public Complex() {
        real = 0.0f;
        imag = 0.0f;
    }

    public Complex(float re, float im) {
        real = re;
        imag = im;
    }

    public float real() {
        return real;
    }

    public float imag() {
        return imag;
    }
}
