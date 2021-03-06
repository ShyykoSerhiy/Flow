package com.github.shyykoserhiy.flow;


/**
 * User: shyyko
 * Date: 02.02.14
 * Time: 19:14
 */
public class Point implements Cloneable {
    private double x;
    private double y;

    public Point(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public Point add(Point point) {
        return new Point(this.getX() + point.getX(), this.getY() + point.getY());
    }

    public Point minus(Point point) {
        return new Point(this.getX() - point.getX(), this.getY() - point.getY());
    }

    public Point divide(double scalar) {
        return new Point(this.getX() / scalar, this.getY() / scalar);
    }

    public Point multiply(double scalar) {
        return new Point(this.getX() * scalar, this.getY() * scalar);
    }

    public double multiply(Point point) {
        return this.getX() * point.getX() + this.getY() * point.getY();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Point point = (Point) o;

        if (Double.compare(point.x, x) != 0) return false;
        if (Double.compare(point.y, y) != 0) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        temp = Double.doubleToLongBits(x);
        result = (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(y);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    @Override
    public String toString() {
        return "Point{" +
                "x=" + x +
                ", y=" + y +
                '}';
    }

    public static boolean intersects(Point firstStartPoint, Point firstEndPoint, Point secondStartPoint, Point secondEndPoint) {
        StraightLine line = new StraightLine(firstEndPoint, secondEndPoint);
        if (!line.projectionIsOnLine(secondStartPoint) && !line.projectionIsOnLine(secondEndPoint)){
            return false;
        }

        Point crossPoint = null;
        Point firstVector = firstEndPoint.minus(firstStartPoint);
        Point secondVector = secondEndPoint.minus(secondStartPoint);

        //computing equations
        double a1 = -firstVector.y;
        double b1 = +firstVector.x;
        double d1 = -(a1 * firstStartPoint.x + b1 * firstStartPoint.y);

        double a2 = -secondVector.y;
        double b2 = +secondVector.x;
        double d2 = -(a2 * secondStartPoint.x + b2 * secondStartPoint.y);

        double seg1Line2Start = a2 * firstStartPoint.x + b2 * firstStartPoint.y + d2;
        double seg1Line2End = a2 * firstEndPoint.x + b2 * firstEndPoint.y + d2;

        double seg2Line1Start = a1 * secondStartPoint.x + b1 * secondStartPoint.y + d1;
        double seg2Line1End = a1 * secondEndPoint.x + b1 * secondEndPoint.y + d1;

        if (seg1Line2Start * seg1Line2End >= 0 || seg2Line1Start * seg2Line1End >= 0)
            return false;

        double u = seg1Line2Start / (seg1Line2Start - seg1Line2End);
        crossPoint = firstStartPoint.add(firstVector.multiply(u)); //this is crossed point

        return true;
    }

    /**
     * Finds mirror point for pointToReflect
     *
     * @param startPoint     beginning of the line
     * @param endPoint       end of the line
     * @param pointToReflect point to reflect against line
     * @return reflected point
     */
    public static Point reflect(Point startPoint, Point endPoint, Point pointToReflect) {
        Point normal = new Point(startPoint.y - endPoint.y, endPoint.x - startPoint.x);
        double length = distance(new Point(0, 0), normal);
        normal = normal.multiply(1 / length);
        StraightLine lLine = new StraightLine(startPoint, endPoint);
        double lDot2 = lLine.distanceToPoint(pointToReflect);
        Point lResult = pointToReflect.minus(normal.multiply(lDot2));
        if (lLine.f(lResult) * lLine.f(pointToReflect) < 0)
            return lResult;
        lResult = pointToReflect.add(normal.multiply(lDot2));
        return lResult;
    }

    public static double distance(double x, double y) {
        return Math.sqrt(x * x + y * y);
    }

    public static double distance(Point A, Point B) {
        return distance(A.x - B.x, A.y - B.y);
    }

    public double Norm() {
        return distance(new Point(0, 0), this);
    }

	@Override
	public Point clone() throws CloneNotSupportedException {
		return new Point(x,y);
	}

	public static class StraightLine {
        Point pA;
        Point pB;

        double attributeA;
        double attributeB;
        double attributeC;

        public StraightLine(Point pA, Point pB) {
            attributeA = pB.y - pA.y;
            attributeB = pA.x - pB.x;
            attributeC = pA.y * pB.x - pA.x * pB.y;
            this.pA = pA;
            this.pB = pB;
        }

        public double f(double pX, double pY) {
            return attributeA * pX + attributeB * pY + attributeC;
        }

        public double f(Point pPoint) {
            return f(pPoint.x, pPoint.y);
        }

        public double distanceToPoint(Point pPoint) {
            if (projectionIsOnLine(pPoint)){
                return rawDistanceToPoint(pPoint);
            } else {
                return 1000000; //todo
            }
        }

        public double rawDistanceToPoint(Point pPoint){
            return Math.abs(f(pPoint)) / Math.sqrt(attributeA * attributeA + attributeB * attributeB);
        }

        public boolean projectionIsOnLine(Point point) {
            double distance = rawDistanceToPoint(point);
            Point normal = new Point(pA.getY() - pB.getY(), pA.getX() - pB.getX());
            normal.divide(distance(normal, new Point(0, 0)));
            Point pointOnLine = point.add(normal.multiply(distance));
            if (Math.abs(f(pointOnLine)) > 0.01) {
                normal = normal.multiply(-1);
                pointOnLine = point.add(normal.multiply(distance));
            }
            if (((pA.getX() <= pointOnLine.getX() && pointOnLine.getX() <= pB.getX()) || ((pA.getX() >= pointOnLine.getX() && pointOnLine.getX() >= pB.getX())) &&
                    ((pA.getY() <= pointOnLine.getY() && pointOnLine.getY() <= pB.getY()) || pA.getY() >= pointOnLine.getY() && pointOnLine.getY() >= pB.getY()))) {
                return true;
            }
            return false;
        }

        public Point pointOnLine() {
            if (attributeA != 0) {
                return new Point(attributeC / attributeA, 0);
            } else if (attributeB != 0) {
                return new Point(0, attributeC / attributeB);
            } else throw new RuntimeException("Invalid line");
        }

        public Point normalToPoint(Point pPoint) {
            Point lPointOnLine = pointOnLine();
            Point lNormal = new Point(attributeA, attributeB);
            lNormal = lNormal.multiply(1 / lNormal.Norm());
            if (distanceToPoint(pPoint) == 0)
                return lNormal;
            if (f(lPointOnLine.add(lNormal)) * f(pPoint) > 0)
                return lNormal;
            else return lNormal.multiply(-1);
        }

        public Point awayFromLine(Point pPoint, double pDistanceAway) {
            double lDistanceToPoint = distanceToPoint(pPoint);
            if (lDistanceToPoint < pDistanceAway) {
                pPoint = pPoint.add((normalToPoint(pPoint)).multiply(pDistanceAway - lDistanceToPoint));
                return pPoint;
            } else return pPoint;
        }
    }
}
