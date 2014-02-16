package shape;

import draw.Point;

import java.util.List;

/**
 * User: shyyko
 * Date: 08.02.14
 * Time: 14:18
 */
public abstract class Shape {
    /**
     * List of points which describe shape
     */
    protected List<Point> listOfPoints;
    /**
     * All points
     */
    protected List<Point> allPoints;
    /**
     * Kolokation points
     */
    protected List<Point> kolokPoints;
    /**
     * Normals in kolocation points
     */
    protected List<Point> normal;
    /**
     * Points from which flying points will take off
     */
    protected List<Point> pointsOfSeparation;

    public List<Point> getListOfPoints() {
        return listOfPoints;
    }

    public List<Point> getAllPoints() {
        return allPoints;
    }

    public List<Point> getKolokPoints() {
        return kolokPoints;
    }

    public List<Point> getNormal() {
        return normal;
    }

    public List<Point> getPointsOfSeparation() {
        return pointsOfSeparation;
    }
}
