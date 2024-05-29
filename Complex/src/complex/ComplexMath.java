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
public final class ComplexMath {

    private ComplexMath() {
    }

    public static Complex add(Complex a, Complex b) {
        return new Complex(a.real() + b.real(), a.imag() + b.imag());
    }

    public static Complex minus(Complex a, Complex b) {
        return new Complex(a.real() - b.real(), a.imag() - b.imag());
    }

    public static Complex times(Complex a, Complex b) {
        return new Complex(a.real() * b.real() - a.imag() * b.imag(),
                a.real() * b.imag() + a.imag() * b.real());
    }

    public static Complex times(Complex a, float alpha) {
        return new Complex(a.real() * alpha, a.imag() * alpha);
    }

    public static Complex timesConjugate(Complex a, Complex b) {
        return new Complex(a.real() * b.real() + a.imag() * b.imag(),
                a.imag() * b.real() - a.real() * b.imag());
    }

    public static Complex conjugate(Complex a) {
        return new Complex(a.real(), -a.imag());
    }

    public static Complex divide(Complex a, float b) {
        return new Complex(a.real() / b, a.imag() / b);
    }

    public Complex divide(Complex a, Complex b) {
        float m = b.real() * b.real() + b.imag() * b.imag();
        return new Complex((a.real() * b.real() + a.imag() * b.imag()) / m,
                (a.imag() * b.real() - a.real() * b.imag()) / m);
    }

    public static Complex cneg(Complex a) {
        return new Complex(-a.real(), -a.imag());
    }

    public static float csqr(Complex a) {
        return (a.real() * a.real()) + (a.imag() * a.imag());
    }

    public static float cabs(Complex a) {
        return (float) Math.sqrt(csqr(a));
    }

    public static Complex cexp(Complex a) {
        if (a.real() == 0.0f) {
            return new Complex((float) Math.cos(a.imag()), (float) Math.sin(a.imag()));
        } else {
            float expf = (float) Math.exp(a.real());

            return new Complex(expf * (float) Math.cos(a.imag()), expf * (float) Math.sin(a.imag()));
        }
    }

    public static float carg(Complex a) {
        return (float) Math.atan2(a.imag(), a.real() + 1E-12f);
    }

    public static Complex normalize(Complex a) {
        float mag = cabs(a);
        return new Complex(a.real() / mag, a.imag() / mag);
    }
}
