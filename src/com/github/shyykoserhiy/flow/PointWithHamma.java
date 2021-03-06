package com.github.shyykoserhiy.flow;

/**
 * User: shyyko
 * Date: 16.02.14
 * Time: 12:42
 */
public class PointWithHamma extends Point {
    private double hamma;

    public PointWithHamma(Point point, double hamma) {
        super(point.getX(), point.getY());
        this.hamma = hamma;
    }

    public PointWithHamma(double x, double y, double hamma) {
        super(x, y);
        this.hamma = hamma;
    }

    public double getHamma() {
        return hamma;
    }

    public void setHamma(double hamma) {
        this.hamma = hamma;
    }

    public boolean equalsByCoordinates(Point point){
        return super.equals(point);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        PointWithHamma that = (PointWithHamma) o;

        if (Double.compare(that.hamma, hamma) != 0) return false;

        return true;
    }

	@Override
	public PointWithHamma clone() throws CloneNotSupportedException {
		return new PointWithHamma(super.clone(), hamma);
	}
}
