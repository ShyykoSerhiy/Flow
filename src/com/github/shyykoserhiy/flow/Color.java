package com.github.shyykoserhiy.flow;

/**
 * User: shyyko
 * Date: 08.02.14
 * Time: 16:47
 */
public class Color {
    private float r;
    private float g;
    private float b;

    public Color(float r, float g, float b) {
        this.r = r;
        this.g = g;
        this.b = b;
    }

    public float getR() {
        return r;
    }

    public void setR(float r) {
        this.r = r;
    }

    public float getG() {
        return g;
    }

    public void setG(float g) {
        this.g = g;
    }

    public float getB() {
        return b;
    }

    public void setB(float b) {
        this.b = b;
    }

    public Color minus(Color color) {
        return new Color(
                getR() - color.getR(),
                getG() - color.getG(),
                getB() - color.getB()
        );
    }

    public Color plus(Color color) {
        return new Color(
                getR() + color.getR(),
                getG() + color.getG(),
                getB() + color.getB()
        );
    }

    public Color divide(float byScalar) {
        return new Color(
                getR() / byScalar,
                getG() / byScalar,
                getB() / byScalar
        );
    }

    public Color multiply(float byScalar) {
        return new Color(
                getR() * byScalar,
                getG() * byScalar,
                getB() * byScalar
        );
    }
}
