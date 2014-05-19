package solver;

import draw.Point;
import draw.PointWithHamma;
import shape.ShShape;
import shape.Shape;

import java.util.*;

import static java.lang.Math.*;
import static java.lang.Math.abs;

/**
 * Created with IntelliJ IDEA.
 * User: shyyko
 * Date: 08.02.14
 * Time: 15:17
 */
public class Solver {
    public static final double DELTA = 0.001;
    private static final int MAXIMUM_OF_FLYING_POINTS_PER_SEPARATION = 2000000;

    private Shape shape;
    /**
     * Is vector which has length 1 i.e. (x^2 +y^2)^(1/2) == 1
     */
    private Point vEternity;
	private double vEternitySquare;

    /**
     * Points
     */
    private Map<PointWithHamma, List<PointWithHamma>> flyingPointsWithHama;
	private List<Pair<PointWithHamma>> dipols = new ArrayList<Pair<PointWithHamma>>();
    private List<Pair<PointWithHamma>> derivativeDipols = new ArrayList<Pair<PointWithHamma>>();
    private double alpha;
	/**
	 * time spent on last move.
	 */
	private double lastMoveSpentTime;
	/**
	 * Contains all states of shapes (cloned ones)
	 */
	private List<Shape> allShapes = new ArrayList<Shape>();

    public Solver(Shape shape) {
        this.shape = shape;
        alpha = 0; //todo
        this.vEternity = new Point(cos(alpha), sin(alpha)); //todo move to constructor parameters
		this.vEternitySquare = vEternity.multiply(vEternity);

        flyingPointsWithHama = new HashMap<PointWithHamma, List<PointWithHamma>>();
        for (PointWithHamma pointWithHamma : shape.getPointsOfSeparation()) {
            flyingPointsWithHama.put(pointWithHamma, new ArrayList<PointWithHamma>());
        }
    }

    public Map<PointWithHamma, List<PointWithHamma>> getFlyingPointsWithHama() {
        return flyingPointsWithHama;
    }

    private boolean firstSolve = true;

    public void solve() {
        // adding flying points (pointsWithHamma)
        if (!firstSolve) {
            for (PointWithHamma pointWithHamma : shape.getPointsOfSeparation()) {
                PointWithHamma flyingPoint = new PointWithHamma(pointWithHamma.getX(), pointWithHamma.getY(), pointWithHamma.getHamma());
                List<PointWithHamma> pointWithHammaList = flyingPointsWithHama.get(pointWithHamma);
                if (pointWithHammaList.size() == MAXIMUM_OF_FLYING_POINTS_PER_SEPARATION) {
                    pointWithHammaList.remove(0);
                }
                pointWithHammaList.add(flyingPoint);
            }
            movePoints();
        } else {
            firstSolve = false;
        }

        List<List<Double>> matrix = new ArrayList<List<Double>>();
        List<Double> rightPartOfSystem = new ArrayList<Double>();
        int amountOfCollocPoints = shape.getKolokPoints().size();

        for (int i = 0; i < amountOfCollocPoints; i++) {
            Point kolokPoint = shape.getKolokPoints().get(i);
            Point normalInKolocPoint = shape.getNormal().get(i);

            List<Double> rowOfMatrix = new ArrayList<Double>(amountOfCollocPoints);
            for (int j = 0; j < amountOfCollocPoints + 1; j++) {
                Point vJ = computeVj(kolokPoint, shape.getAllPoints().get(j));
                double coefficientToAdd = vJ.multiply(normalInKolocPoint);
                rowOfMatrix.add(coefficientToAdd);
            }
            matrix.add(rowOfMatrix);

            double rightPart = -vEternity.multiply(shape.getNormal().get(i));
            for (List<PointWithHamma> list : flyingPointsWithHama.values()) {
                for (PointWithHamma pointWithHamma : list) {
                    rightPart -= pointWithHamma.getHamma() * computeVj(kolokPoint, pointWithHamma).multiply(normalInKolocPoint);
                }
            }

            rightPartOfSystem.add(rightPart);
        }

        matrix.add(new ArrayList<Double>(Collections.nCopies(amountOfCollocPoints + 1, 1.0)));
        double lastRightPart = 0;
        for (List<PointWithHamma> list : flyingPointsWithHama.values()) {
            for (PointWithHamma pointWithHamma : list) {
                lastRightPart -= pointWithHamma.getHamma();
            }
        }


        rightPartOfSystem.add(lastRightPart);

        List<Double> hamma = gaussMethod(matrix, rightPartOfSystem);
        for (int i = 0; i < hamma.size(); i++) {
            shape.getAllPoints().get(i).setHamma(hamma.get(i));
        }

        int amountBySide = (shape.getAllPoints().size() - 6) / 4;
        shape.getAllPoints().get(shape.getAllPoints().size() - 1)
                .setHamma(shape.getAllPoints().get(2 + amountBySide + (amountBySide / 2) + 1).getHamma());

		try {
			allShapes.add(shape.clone());
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		//splitting to dipols
		splitContourToDipols();
        splitContourToDerivativeDipols();
    }

	private void splitContourToDipols(){
		List<List<PointWithHamma>> splitContour = new ArrayList<List<PointWithHamma>>();

		splitContour.add(getSplitFlyingPointsPart(0));     //////////////////////////////
		splitContour.add(getSplitFlyingPointsPart(1));      //////////////////////////////
		splitContour.add(getSplitFlyingPointsPart(2));      //////////////////////////////
		splitContour.add(getSplitFlyingPointsPart(3));      //////////////////////////////
		splitContour.add(getSplitFlyingPointsPart(4));     //////////////////////////////
		PointWithHamma connectingPoint0 = splitContour.get(0).get(splitContour.get(0).size() - 1);
		PointWithHamma connectingPoint1 = splitContour.get(1).get(splitContour.get(1).size() - 1);
		PointWithHamma connectingPoint2 = splitContour.get(2).get(splitContour.get(2).size() - 1);
		PointWithHamma connectingPoint3 = splitContour.get(3).get(splitContour.get(3).size() - 1);
		PointWithHamma connectingPoint4 = splitContour.get(4).get(splitContour.get(4).size() - 1);

		PointWithHamma firstPoint = shape.getPointsOfSeparation().get(4);
		List<PointWithHamma> splitPart = new ArrayList<PointWithHamma>();
		splitPart.add(new PointWithHamma(firstPoint, -connectingPoint4.getHamma() + firstPoint.getHamma()));
		double sumHamma = splitPart.get(0).getHamma();
		for (int i = shape.getAllPoints().indexOf(shape.getPointsOfSeparation().get(4)) + 1;
			 i < shape.getAllPoints().size(); i++) {
			PointWithHamma pointWithHamma = shape.getAllPoints().get(i);
			PointWithHamma pointToAdd = new PointWithHamma(pointWithHamma, pointWithHamma.getHamma());
			splitPart.add(pointToAdd);
			sumHamma += pointToAdd.getHamma();
		}
		PointWithHamma lastPoint = splitPart.get(splitPart.size() - 1);
		sumHamma -= lastPoint.getHamma();
		lastPoint.setHamma(-sumHamma);
		splitContour.add(splitPart);    //////////////////////////////
		PointWithHamma connectingPoint5 = splitContour.get(5).get(splitContour.get(5).size() - 1);

		splitPart = new ArrayList<PointWithHamma>();
		for (int i = 0; i <= shape.getAllPoints().indexOf(shape.getPointsOfSeparation().get(3)); i++) {
			PointWithHamma pointWithHamma = shape.getAllPoints().get(i);
			if (pointWithHamma.equalsByCoordinates(connectingPoint0)) {
				splitPart.add(new PointWithHamma(pointWithHamma, pointWithHamma.getHamma() - connectingPoint0.getHamma()));
			} else if (pointWithHamma.equalsByCoordinates(connectingPoint1)) {
				splitPart.add(new PointWithHamma(pointWithHamma, pointWithHamma.getHamma() - connectingPoint1.getHamma()));
			} else if (pointWithHamma.equalsByCoordinates(connectingPoint2)) {
				splitPart.add(new PointWithHamma(pointWithHamma, pointWithHamma.getHamma() - connectingPoint2.getHamma()));
			} else if (pointWithHamma.equalsByCoordinates(connectingPoint3)) {
				splitPart.add(new PointWithHamma(pointWithHamma, pointWithHamma.getHamma() - connectingPoint3.getHamma()));
			} else if (pointWithHamma.equalsByCoordinates(connectingPoint5)) {
				splitPart.add(new PointWithHamma(pointWithHamma, pointWithHamma.getHamma() - connectingPoint5.getHamma()));
			} else {
				splitPart.add(new PointWithHamma(pointWithHamma, pointWithHamma.getHamma()));
			}
		}
		splitContour.add(splitPart);

		double summ = 0;
		for (List<PointWithHamma> pointWithHammaList : flyingPointsWithHama.values()) {
			for (PointWithHamma pointWithHamma : pointWithHammaList) {
				summ += pointWithHamma.getHamma();
			}
		}

		for (int i = 0; i < splitContour.size(); i++) {
			List<PointWithHamma> list = splitContour.get(i);
			double sum = 0;
			for (PointWithHamma pointWithHamma : list) {
				sum += pointWithHamma.getHamma();
			}
			if (abs(sum) > 0.01 && (abs(summ - sum) > 0.01)) {
				throw new RuntimeException("invalid sum");
			}
		}


		dipols = new ArrayList<Pair<PointWithHamma>>();
		for (List<PointWithHamma> pointWithHammaList : splitContour) {
			double oldHama = 0;
			for (int i = 0; i < pointWithHammaList.size() - 1; i++) {
				PointWithHamma firstPointWithHamma = pointWithHammaList.get(i);
				PointWithHamma secondPointWithHamma = pointWithHammaList.get(i + 1);

				PointWithHamma newFirstPointWithHamma = new PointWithHamma(firstPointWithHamma,
						firstPointWithHamma.getHamma() + oldHama);
				PointWithHamma newSecondPointWithHamma = new PointWithHamma(secondPointWithHamma, -newFirstPointWithHamma.getHamma());
				oldHama = newFirstPointWithHamma.getHamma();
				Pair<PointWithHamma> dipol = new Pair<PointWithHamma>(newFirstPointWithHamma, newSecondPointWithHamma);
				dipols.add(dipol);
			}
		}
	}

    private void splitContourToDerivativeDipols(){
        if (allShapes.size() < 2){
            return;
        }
        Shape currentShape = allShapes.get(allShapes.size()-1);
        Shape previousShape = allShapes.get(allShapes.size()-2);

        List<List<PointWithHamma>> splitContour = new ArrayList<List<PointWithHamma>>();

        splitContour.add(getSplitDerivativeFlyingPointsPart(0));     //////////////////////////////
        splitContour.add(getSplitDerivativeFlyingPointsPart(1));      //////////////////////////////
        splitContour.add(getSplitDerivativeFlyingPointsPart(2));      //////////////////////////////
        splitContour.add(getSplitDerivativeFlyingPointsPart(3));      //////////////////////////////
        splitContour.add(getSplitDerivativeFlyingPointsPart(4));     //////////////////////////////
        PointWithHamma connectingPoint0 = splitContour.get(0).get(splitContour.get(0).size() - 1);
        PointWithHamma connectingPoint1 = splitContour.get(1).get(splitContour.get(1).size() - 1);
        PointWithHamma connectingPoint2 = splitContour.get(2).get(splitContour.get(2).size() - 1);
        PointWithHamma connectingPoint3 = splitContour.get(3).get(splitContour.get(3).size() - 1);
        PointWithHamma connectingPoint4 = splitContour.get(4).get(splitContour.get(4).size() - 1);

        PointWithHamma firstPoint = currentShape.getPointsOfSeparation().get(4);
        PointWithHamma oldFirstPoint = previousShape.getPointsOfSeparation().get(4);
        List<PointWithHamma> splitPart = new ArrayList<PointWithHamma>();
        splitPart.add(new PointWithHamma(firstPoint, -connectingPoint4.getHamma() + ((firstPoint.getHamma()- oldFirstPoint.getHamma())/lastMoveSpentTime)));
        double sumHamma = splitPart.get(0).getHamma();
        for (int i = currentShape.getAllPoints().indexOf(currentShape.getPointsOfSeparation().get(4)) + 1;
             i < currentShape.getAllPoints().size(); i++) {
            PointWithHamma pointWithHamma = currentShape.getAllPoints().get(i);
            PointWithHamma oldPointWithHammma = previousShape.getAllPoints().get(i);
            double derivativeHamma = (pointWithHamma.getHamma() - oldPointWithHammma.getHamma()) / lastMoveSpentTime;
            PointWithHamma pointToAdd = new PointWithHamma(pointWithHamma, derivativeHamma);
            splitPart.add(pointToAdd);
            sumHamma += derivativeHamma;
        }
        PointWithHamma lastPoint = splitPart.get(splitPart.size() - 1);
        sumHamma -= lastPoint.getHamma();
        lastPoint.setHamma(-sumHamma);
        splitContour.add(splitPart);    //////////////////////////////
        PointWithHamma connectingPoint5 = splitContour.get(5).get(splitContour.get(5).size() - 1);

        splitPart = new ArrayList<PointWithHamma>();
        for (int i = 0; i <= currentShape.getAllPoints().indexOf(currentShape.getPointsOfSeparation().get(3)); i++) {
            PointWithHamma pointWithHamma = currentShape.getAllPoints().get(i);
            PointWithHamma oldPointWithHammma = previousShape.getAllPoints().get(i);
            double derivativeHamma = (pointWithHamma.getHamma() - oldPointWithHammma.getHamma())/lastMoveSpentTime;

            if (pointWithHamma.equalsByCoordinates(connectingPoint0)) {
                splitPart.add(new PointWithHamma(pointWithHamma, derivativeHamma - connectingPoint0.getHamma()));
            } else if (pointWithHamma.equalsByCoordinates(connectingPoint1)) {
                splitPart.add(new PointWithHamma(pointWithHamma, derivativeHamma - connectingPoint1.getHamma()));
            } else if (pointWithHamma.equalsByCoordinates(connectingPoint2)) {
                splitPart.add(new PointWithHamma(pointWithHamma, derivativeHamma - connectingPoint2.getHamma()));
            } else if (pointWithHamma.equalsByCoordinates(connectingPoint3)) {
                splitPart.add(new PointWithHamma(pointWithHamma, derivativeHamma - connectingPoint3.getHamma()));
            } else if (pointWithHamma.equalsByCoordinates(connectingPoint5)) {
                splitPart.add(new PointWithHamma(pointWithHamma, derivativeHamma - connectingPoint5.getHamma()));
            } else {
                splitPart.add(new PointWithHamma(pointWithHamma, derivativeHamma));
            }
        }
        splitContour.add(splitPart);

       /* double summ = 0; //todo
        for (List<PointWithHamma> pointWithHammaList : flyingPointsWithHama.values()) {
            for (PointWithHamma pointWithHamma : pointWithHammaList) {
                summ += pointWithHamma.getHamma();
            }
        }

        for (int i = 0; i < splitContour.size(); i++) {
            List<PointWithHamma> list = splitContour.get(i);
            double sum = 0;
            for (PointWithHamma pointWithHamma : list) {
                sum += pointWithHamma.getHamma();
            }
            if (abs(sum) > 0.01 && (abs(summ - sum) > 0.01)) {
                throw new RuntimeException("invalid sum");
            }
        }*/


        derivativeDipols = new ArrayList<Pair<PointWithHamma>>();
        for (List<PointWithHamma> pointWithHammaList : splitContour) {
            double oldHama = 0;
            for (int i = 0; i < pointWithHammaList.size() - 1; i++) {
                PointWithHamma firstPointWithHamma = pointWithHammaList.get(i);
                PointWithHamma secondPointWithHamma = pointWithHammaList.get(i + 1);

                PointWithHamma newFirstPointWithHamma = new PointWithHamma(firstPointWithHamma,
                        firstPointWithHamma.getHamma() + oldHama);
                PointWithHamma newSecondPointWithHamma = new PointWithHamma(secondPointWithHamma, -newFirstPointWithHamma.getHamma());
                oldHama = newFirstPointWithHamma.getHamma();
                Pair<PointWithHamma> dipol = new Pair<PointWithHamma>(newFirstPointWithHamma, newSecondPointWithHamma);
                derivativeDipols.add(dipol);
            }
        }
    }

	private List<PointWithHamma> getSplitFlyingPointsPart(int numberOfSeparationPoint) {
		PointWithHamma pointOfSeparation = shape.getPointsOfSeparation().get(numberOfSeparationPoint);
		List<PointWithHamma> splitPart = new ArrayList<PointWithHamma>();
		double sumHama = 0;
		for (PointWithHamma pointWithHamma : flyingPointsWithHama.get(pointOfSeparation)) {
			splitPart.add(new PointWithHamma(pointWithHamma, pointWithHamma.getHamma()));
			sumHama += pointWithHamma.getHamma();
		}
		splitPart.add(new PointWithHamma(pointOfSeparation, -sumHama));
		return splitPart;
	}

    private List<PointWithHamma> getSplitDerivativeFlyingPointsPart(int numberOfSeparationPoint) {
        PointWithHamma pointOfSeparation = shape.getPointsOfSeparation().get(numberOfSeparationPoint);
        List<PointWithHamma> splitPart = new ArrayList<PointWithHamma>();
        double sumHama = 0;
        List<PointWithHamma> flyingPoints = flyingPointsWithHama.get(pointOfSeparation);
        PointWithHamma lastSeparated = flyingPoints.get(flyingPoints.size()-1);
        splitPart.add(new PointWithHamma(lastSeparated, lastSeparated.getHamma()/lastMoveSpentTime));
        sumHama += lastSeparated.getHamma()/lastMoveSpentTime;

        splitPart.add(new PointWithHamma(pointOfSeparation, -sumHama));
        return splitPart;
    }

    /**
     * Moving flying points with intersection check
     */
    public void movePoints() {
        if (flyingPointsWithHama.get(shape.getPointsOfSeparation().get(0)).isEmpty()) {
            return;
        }

        double maxSpeed = 0;
        for (List<PointWithHamma> flyingPoints : flyingPointsWithHama.values()) {
            PointWithHamma pointWithHamma = flyingPoints.get(flyingPoints.size() - 1);
            Point speed = computeV(pointWithHamma);
            double speedDouble = Point.distance(speed, new Point(0, 0));
            if (maxSpeed < speedDouble) {
                maxSpeed = speedDouble;
            }
        }
		lastMoveSpentTime = shape.getDelta() / maxSpeed;

        for (int j = 0; j < shape.getPointsOfSeparation().size(); j++) {
            PointWithHamma separationPoint = shape.getPointsOfSeparation().get(j);
            List<PointWithHamma> flyingPoints = flyingPointsWithHama.get(separationPoint);
            int[] pointsToSkip;
            if (j == 1 || j == 2) {
                pointsToSkip = new int[]{3, 4};
            } else {
                pointsToSkip = new int[]{3};
            }
            for (int i = 0; i < flyingPoints.size(); i++) {
                PointWithHamma oldPoint = flyingPoints.get(i);
                Point newPointCoordinate = null;
                newPointCoordinate = computeV(oldPoint).multiply(lastMoveSpentTime).add(oldPoint);
                /*//todo live hack
                if (Point.distance(oldPoint, newPointCoordinate) > shape.getDelta()*1.5){
                   newPointCoordinate = computeV(oldPoint).multiply(shape.getDelta()/Point.distance(oldPoint, newPointCoordinate)).add(oldPoint);
                }
                //todo live hack*/
                if (j==2 && i == 0){
                    System.out.println(newPointCoordinate);
                }
                if (i < flyingPoints.size() - 1) { // there is no need for additional checks for new flying point
                    newPointCoordinate = reflectIsCrossCarcases(oldPoint, newPointCoordinate, pointsToSkip);
                    if (j==2 && i == 0){
                        System.out.println(newPointCoordinate);
                    }
                    newPointCoordinate = awayFromSegmentationPoint(newPointCoordinate, pointsToSkip);
                    if (j==2 && i == 0){
                        System.out.println(newPointCoordinate);
                    }
                }


                flyingPoints.set(i, new PointWithHamma(newPointCoordinate.getX(), newPointCoordinate.getY(),
                        oldPoint.getHamma()));
            }

        }
    }

    public double getPhi(Point point) {
        double toReturn = point.getX() * cos(alpha) + point.getY() * sin(alpha);

        double sum = 0;
        for (int i = 0; i < shape.getAllPoints().size(); i++) {
            sum += shape.getAllPoints().get(i).getHamma() * atan2(
                    point.getY() - shape.getAllPoints().get(i).getY(),
                    point.getX() - shape.getAllPoints().get(i).getX()
            );
        }

        for (List<PointWithHamma> pointWithHammaList : flyingPointsWithHama.values()) {
            for (PointWithHamma pointWithHamma : pointWithHammaList) {
                sum += pointWithHamma.getHamma() * atan2(
                        point.getY() - pointWithHamma.getY(),
                        point.getX() - pointWithHamma.getX()
                );
            }
        }
        toReturn += 1 / (2 * PI) * sum;
        return toReturn;
    }

    //todo this method depends on shShape and is a piece of %#&@*!. Needs to be refactored
    public double getPhiTrueDipol(Point point) {
        double phi = point.getX() * cos(alpha) + point.getY() * sin(alpha);
        double x = point.getX();
        double y = point.getY();

        for (Pair<PointWithHamma> dipol : dipols) {
            double x1 = dipol.getFirstValue().getX();
            double y1 = dipol.getFirstValue().getY();
            double x2 = dipol.getSecondValue().getX();
            double y2 = dipol.getSecondValue().getY();

            phi += (dipol.getFirstValue().getHamma() / (2 * PI)) *
                    ((y2 - y1) * (x - (x2 + x1) / 2) - (x2 - x1) * (y - (y2 + y1) / 2)) /
                    (pow(x - (x2 + x1) / 2, 2) + pow(y - (y2 + y1) / 2, 2));
        }

        /*splitContour.add(splitPart);

        splitPart = new ArrayList<PointWithHamma>();
        splitPart.add(new PointWithHamma(pointOfSeparation, -sumHama + pointOfSeparation.getHamma()));
*/
        return phi;
    }

	public double getCp(Point point) {
		Point vInPoint = computeV(point);
		double cp = 1 - vInPoint.multiply(vInPoint)/vEternitySquare - (2/vEternitySquare) * getPhiDerivative(point);
		return cp;
	}

	@SuppressWarnings("unchecked")
	private double getPhiDerivative(Point point) {
        if (allShapes.size() < 2){
            return 0.0;
        }
		double phiDerivative = 0;

        double phiDerivativeDipol = getPhiTrueDipolDerivative(point);

        double phiDerivativeConv = 0;

        for (List<PointWithHamma> points : flyingPointsWithHama.values()){
            for (int i = 0; i < points.size(); i++) {
                PointWithHamma pointWithHamma = points.get(i);
                phiDerivativeConv += pointWithHamma.getHamma() * computeVj(point, pointWithHamma).multiply(computeV(pointWithHamma));
            }
        }

        phiDerivative = phiDerivativeDipol - phiDerivativeConv;

        return phiDerivative;


		/*if (allShapes.size() > 1) { //from second step it is possible to calculate phi derivative
			//todo move this to solve method for performance
			List<PointWithHamma> derivativeContourPoints = new ArrayList<PointWithHamma>();

			Shape currentShape = allShapes.get(allShapes.size()-1);
			Shape previousShape = allShapes.get(allShapes.size()-2);

			double hammaDerivativeSum = 0;
			for (int i = 0; i < currentShape.getAllPoints().size()-1; i++) {
				PointWithHamma currentPoint = currentShape.getAllPoints().get(i);
				PointWithHamma previousPoint = previousShape.getAllPoints().get(i);

				double hammaDerivative = (currentPoint.getHamma() - previousPoint.getHamma()) / lastMoveSpentTime;
				hammaDerivativeSum += hammaDerivative;

				derivativeContourPoints.add(new PointWithHamma(currentPoint, hammaDerivativeSum));
			}
			derivativeContourPoints.add(new PointWithHamma(
					currentShape.getAllPoints().get(currentShape.getAllPoints().size() - 1), -hammaDerivativeSum
			));

			//todo refactor this whole "clone" thing
			Map<PointWithHamma, List<PointWithHamma>> derivativeFlyingPoints = cloneMap(flyingPointsWithHama);
			for (List<PointWithHamma> points : derivativeFlyingPoints.values()) {
				PointWithHamma lastSeparatedPoint = points.get(points.size() - 1);
				lastSeparatedPoint.setHamma(lastSeparatedPoint.getHamma() / lastMoveSpentTime);
			}

			for (int i = 0; i < currentShape.getAllPoints().size()-1; i++) {
				PointWithHamma pointWithHamma = currentShape.getAllPoints().get(i);
				PointWithHamma nextPointWithHamma = currentShape.getAllPoints().get(i+1);
				double xWithCap = 0.5*(nextPointWithHamma.getX() + pointWithHamma.getX());
				double yWithCap = 0.5*(nextPointWithHamma.getY() + pointWithHamma.getY());


				phiDerivative += (((yWithCap - point.getY()) * (nextPointWithHamma.getX() - pointWithHamma.getX()) +
						(point.getX() - xWithCap) * (nextPointWithHamma.getY() - pointWithHamma.getY())) /
						(Math.pow(point.getX() - xWithCap, 2) + Math.pow(point.getY() - yWithCap, 2))) *
						(pointWithHamma.getHamma() / (2 * PI));
			}

			for (PointWithHamma keyPoint : derivativeFlyingPoints.keySet()) { //todo ?? is it correct(formulas are not intuitive)
				List<PointWithHamma> list = derivativeFlyingPoints.get(keyPoint);
				PointWithHamma pointWithHamma = list.get(list.size()-1);
				double xWithCap = 0.5 * (keyPoint.getX() + pointWithHamma.getX());
				double yWithCap = 0.5 * (keyPoint.getY() + pointWithHamma.getY());


				phiDerivative += (((yWithCap - point.getY()) * (keyPoint.getX() - pointWithHamma.getX()) +
						(point.getX() - xWithCap) * (keyPoint.getY() - pointWithHamma.getY())) /
						(Math.pow(point.getX() - xWithCap, 2) + Math.pow(point.getY() - yWithCap, 2))) *
						(pointWithHamma.getHamma() / (2 * PI));
			}

			for (List<PointWithHamma> points : derivativeFlyingPoints.values()){
				for (int i = 0; i < points.size()-1; i++) {
					PointWithHamma pointWithHamma = points.get(i);
					phiDerivative -= pointWithHamma.getHamma() * computeVj(point, pointWithHamma).multiply(computeV(pointWithHamma));
				}
			}
		}           */
	}

    public double getPhiTrueDipolDerivative(Point point) {
        double phi = 0;
        double x = point.getX();
        double y = point.getY();

        for (Pair<PointWithHamma> dipol : derivativeDipols) {
            double x1 = dipol.getFirstValue().getX();
            double y1 = dipol.getFirstValue().getY();
            double x2 = dipol.getSecondValue().getX();
            double y2 = dipol.getSecondValue().getY();

            phi += (dipol.getFirstValue().getHamma() / (2 * PI)) *
                    ((y2 - y1) * (x - (x2 + x1) / 2) - (x2 - x1) * (y - (y2 + y1) / 2)) /
                    (pow(x - (x2 + x1) / 2, 2) + pow(y - (y2 + y1) / 2, 2));
        }

        /*splitContour.add(splitPart);

        splitPart = new ArrayList<PointWithHamma>();
        splitPart.add(new PointWithHamma(pointOfSeparation, -sumHama + pointOfSeparation.getHamma()));
*/
        return phi;
    }

	private Map<PointWithHamma, List<PointWithHamma>> cloneMap(Map<PointWithHamma, List<PointWithHamma>> flyingPointsWithHama) {
		Map<PointWithHamma, List<PointWithHamma>> toReturn = new HashMap<PointWithHamma, List<PointWithHamma>>();
		for (PointWithHamma pointWithHamma : flyingPointsWithHama.keySet()){
			toReturn.put(pointWithHamma, ShShape.cloneList(flyingPointsWithHama.get(pointWithHamma)));
		}
		return toReturn;
	}

	public Pair<Double> computeForFlyingPoints(double partialHama, Point point, int firstFlyingPointNumber,
                                               int numberOfLastPoint) {
        Point firstStrangeValue = new Point(0, 0);
        Point secondStrangeValue = new Point(0, 0);
        double toReturn = 0;
        List<PointWithHamma> flyingPoints = flyingPointsWithHama.get(shape.getPointsOfSeparation().get(firstFlyingPointNumber));

        for (int i = firstFlyingPointNumber; i < flyingPoints.size(); i++) {
            PointWithHamma flyingPoint = flyingPoints.get(i);
            partialHama += flyingPoint.getHamma();

            if (i + 1 < flyingPoints.size()) {
                Point neighborPoint = flyingPoints.get(i + 1);
                firstStrangeValue = neighborPoint.minus(flyingPoint);
            } else {
                Point neighborPoint = shape.getListOfPoints().get(numberOfLastPoint);
                firstStrangeValue = neighborPoint.minus(flyingPoint);
            }

            secondStrangeValue = point.minus(flyingPoint);

            double strangeR = Math.pow(secondStrangeValue.getX(), 2) + Math.pow(secondStrangeValue.getY(), 2);
            strangeR = strangeR < DELTA ? DELTA : strangeR;

            double strangeProd = (firstStrangeValue.getY() * secondStrangeValue.getX() - firstStrangeValue.getX() * secondStrangeValue.getY())
                    / (strangeR * (2 * PI));
            strangeProd *= partialHama;

            toReturn += strangeProd;
        }

        return new Pair<Double>(toReturn, partialHama);
    }

    public Pair<Double> computeForContourPoints(double partialHama, Point point,
                                                int firstPointNumber, int lastPointNumber) {
        Point firstStrangeValue = new Point(0, 0);
        Point secondStrangeValue = new Point(0, 0);
        double toReturn = 0;

        for (int i = firstPointNumber; i < lastPointNumber; i += 1) {
            Point pointOnContour = shape.getAllPoints().get(i);
            partialHama += shape.getAllPoints().get(i).getHamma();

            Point neighborPoint = shape.getAllPoints().get(i + 1);
            firstStrangeValue = neighborPoint.minus(pointOnContour);

            secondStrangeValue = point.minus(pointOnContour);

            double strangeR = Math.pow(secondStrangeValue.getX(), 2) + Math.pow(secondStrangeValue.getY(), 2);
            strangeR = strangeR < DELTA ? DELTA : strangeR;

            double strangeProd = (firstStrangeValue.getY() * secondStrangeValue.getX() - firstStrangeValue.getX() * secondStrangeValue.getY())
                    / (strangeR * (2 * PI));
            strangeProd *= partialHama;

            toReturn += strangeProd;
        }

        return new Pair<Double>(toReturn, partialHama);
    }



    public double getPsi(Point point) {
        double toReturn = cos(alpha) * point.getY() - sin(alpha) * point.getX();

        double sum = 0;
        for (int i = 0; i < shape.getAllPoints().size(); i++) {
            sum += shape.getAllPoints().get(i).getHamma() * log(computeRj(point, shape.getAllPoints().get(i)));
        }
        for (List<PointWithHamma> flyingPoints : flyingPointsWithHama.values()) {
            for (PointWithHamma flyingPoint : flyingPoints) {
                sum += flyingPoint.getHamma() * log(computeRj(point, flyingPoint));
            }
        }
        toReturn -= 1 / (2 * PI) * sum;
        return toReturn;
    }

    private Point computeVj(Point firstPoint, Point secondPoint) {
        return new Point(
                getU(firstPoint, secondPoint),
                getV(firstPoint, secondPoint)
        );
    }

    public Point computeV(Point point) {
        Point sum = vEternity.add(new Point(0, 0));
        for (int i = 0; i < shape.getAllPoints().size(); i++) {
            sum = sum.add(computeVj(point, shape.getAllPoints().get(i)).multiply(shape.getAllPoints().get(i).getHamma()));
        }
        for (List<PointWithHamma> flyingPoints : flyingPointsWithHama.values()) {
            for (PointWithHamma flyingPoint : flyingPoints) {
                sum = sum.add(computeVj(point, flyingPoint).multiply(flyingPoint.getHamma()));
            }
        }
        return sum;
    }

    private double getU(Point firstPoint, Point secondPoint) {
        return 1 / (2 * PI) * ((secondPoint.getY() - firstPoint.getY()) / pow(computeRj(firstPoint, secondPoint), 2));
    }

    private double getV(Point firstPoint, Point secondPoint) {
        return 1 / (2 * PI) * ((firstPoint.getX() - secondPoint.getX()) / pow(computeRj(firstPoint, secondPoint), 2));
    }

    private double computeRj(Point firstPoint, Point secondPoint) {
        return max(getSigma(), sqrt(
                pow(firstPoint.getX() - secondPoint.getX(), 2) + pow(firstPoint.getY() - secondPoint.getY(), 2))
        );
    }

    public double getSigma() { //todo implement
        return 0.001;
    }


    private static List<Double> gaussElimination(List<List<Double>> matrix, List<Double> rightSide) {
        // Gaussian elimination with partial pivoting
        int size = rightSide.size();

        for (int p = 0; p < size; p++) {
            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < size; i++) {
                if (abs(matrix.get(i).get(p)) > abs(matrix.get(max).get(p))) {
                    max = i;
                }
            }
            List<Double> temp = matrix.get(p);
            matrix.set(p, matrix.get(max));
            matrix.set(max, temp);
            double t = rightSide.get(p);
            rightSide.set(p, rightSide.get(max));
            rightSide.set(max, t);

            // singular or nearly singular
            if (abs(matrix.get(p).get(p)) <= 1e-10) {
                throw new RuntimeException("Matrix is singular or nearly singular");
            }

            // pivot within A and b
            for (int i = p + 1; i < size; i++) {
                double alpha = matrix.get(i).get(p) / matrix.get(p).get(p);
                rightSide.set(i, rightSide.get(i) - alpha * rightSide.get(p));
                for (int j = p; j < size; j++) {
                    matrix.get(i).set(j, matrix.get(i).get(j) - alpha * matrix.get(p).get(j));
                }
            }
        }

        // back substitution
        Double[] x = new Double[size];
        for (int i = size - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < size; j++) {
                sum += matrix.get(i).get(j) * x[j];
            }
            x[i] = (rightSide.get(i) - sum) / matrix.get(i).get(i);
        }
        return new ArrayList<Double>(Arrays.asList(x));
    }

    public static final float EPS = 0.001f;


    private static List<Double> gaussMethod(List<List<Double>> matrix, List<Double> rightSide) {
        float[][] matrixA = new float[rightSide.size()][rightSide.size()];
        float[] b = new float[rightSide.size()];

        for (int i = 0; i < rightSide.size(); i++) {
            for (int j = 0; j < rightSide.size(); j++) {
                matrixA[i][j] = matrix.get(i).get(j).floatValue();
            }
            b[i] = rightSide.get(i).floatValue();
        }
        gaussMethod(matrixA, b);
        List<Double> toReturn = new ArrayList<Double>(rightSide.size());
        for (float value : b) {
            toReturn.add(Double.valueOf(value));
        }
        return toReturn;
    }

    public static void gaussMethod(float[][] a, float[] b) {
        for (int i = 0; i < b.length - 1; i++) {
            int l = 0;
            while (Math.abs(a[i][i]) < EPS) {
                l++;
                if (i < b.length - l) {
                    for (int j = 0; j < b.length; j++) {
                        a[i][j] += a[i + l][j];
                    }
                    b[i] += b[i + l];
                } else {
                    for (int j = 0; j < b.length; j++) {
                        a[i][j] += a[i - l][j];
                    }
                    b[i] += b[i - l];
                }
            }

            for (int j = i + 1; j < b.length; j++) {
                if (Math.abs(a[j][i]) > EPS) {
                    final float koef = -a[j][i] / a[i][i];
                    for (int k = 0; k < b.length; k++) {
                        a[j][k] += koef * a[i][k];
                    }
                    b[j] += b[i] * koef;
                }
            }
        }
        for (int j = b.length - 1; j >= 1; j--) {
            for (int i = 0; i < j; i++) {
                if (Math.abs(a[i][j]) > EPS) {
                    final float koef = -a[i][j] / a[j][j];
                    a[i][j] += koef * a[j][j];
                    b[i] += b[j] * koef;
                }
            }
            b[j] /= a[j][j];
            a[j][j] = 1;
        }
        b[0] /= a[0][0];
        a[0][0] = 1;
    }

    /**
     * Reflects new point if it crosses the carcases
     *
     * @param oldPoint
     * @param newPoint
     * @return
     */
    private Point reflectIsCrossCarcases(Point oldPoint, Point newPoint, int[] pointsToSkip) {
        boolean lResult = true;
        int j = 0;
        mainCycle:
        while (j < shape.getListOfPoints().size() - 1) {
            for (int pointToSkip : pointsToSkip) {
                if (j == pointToSkip) {
                    j++;
                    continue mainCycle;
                }
            }
            Point startPoint = shape.getListOfPoints().get(j);
            Point endPoint = shape.getListOfPoints().get(j + 1);
            if (Point.intersects(newPoint, oldPoint, startPoint, endPoint)) {

                return away(startPoint, endPoint, oldPoint, newPoint);
                //newPoint = Point.reflect(shape.getAllPoints().get(j), shape.getAllPoints().get(j + 1), newPoint);
                /*j = 0;
                lResult = false;*/
            } else j++;
        }
        return newPoint;
    }

    public Point away(Point startPoint, Point endPoint, Point oldPoint, Point newPoint) {
        Point normal = new Point(startPoint.getY() - endPoint.getY(), endPoint.getX() - startPoint.getX());
        double distance = Point.distance(normal, new Point(0, 0));
        normal = normal.divide(distance);

        Point.StraightLine line = new Point.StraightLine(startPoint, endPoint);
        double distanceToLine = line.distanceToPoint(newPoint);
        if (oldPoint.equals(newPoint)) {
            distanceToLine = 0;
        }
        Point newPossiblePoint = newPoint.add(normal.multiply(shape.getDelta() + distanceToLine));
        if (!Point.intersects(newPossiblePoint, oldPoint, startPoint, endPoint)) {
            return newPossiblePoint;
        } else {
            normal = normal.multiply(-1);
            return newPoint.add(normal.multiply(shape.getDelta() + distanceToLine));
        }
    }

    public Point awaySeg(Point startPoint, Point endPoint, Point oldPoint, Point newPoint) {
        Point normal = new Point(startPoint.getY() - endPoint.getY(), endPoint.getX() - startPoint.getX());
        double distance = Point.distance(normal, new Point(0, 0));
        normal = normal.divide(distance);

        Point.StraightLine line = new Point.StraightLine(startPoint, endPoint);
        double distanceToLine = line.distanceToPoint(newPoint);
        Point newPossiblePoint = newPoint.add(normal.multiply(shape.getDelta() - distanceToLine));
        if (abs(line.distanceToPoint(newPossiblePoint) - shape.getDelta()) < 0.01 && !Point.intersects(startPoint, endPoint, oldPoint, newPossiblePoint)) {
            return newPossiblePoint;
        } else {
            normal = normal.multiply(-1);
            return newPoint.add(normal.multiply(shape.getDelta() - distanceToLine));
        }
    }

    public static void main(String[] args) {
        Point firstPoint = new Point(2, 2);
        Point secondPoint = new Point(1, 2);
        Point newPoint = new Point(1.2f, 2.04f);
        awaySega(firstPoint, secondPoint, newPoint, newPoint);

    }

    public static Point awaySega(Point startPoint, Point endPoint, Point oldPoint, Point newPoint) {
        double delta = 0.083333;
        Point normal = new Point(startPoint.getY() - endPoint.getY(), endPoint.getX() - startPoint.getX());
        double distance = Point.distance(normal, new Point(0, 0));
        normal = normal.divide(distance);

        Point.StraightLine line = new Point.StraightLine(startPoint, endPoint);
        double distanceToLine = line.distanceToPoint(newPoint);
        Point newPossiblePoint = newPoint.add(normal.multiply(delta - distanceToLine));
        if (abs(line.distanceToPoint(newPossiblePoint) - delta) < 0.01 && !Point.intersects(startPoint, endPoint, oldPoint, newPossiblePoint)) {
            return newPossiblePoint;
        } else {
            normal = normal.multiply(-1);
            return newPoint.add(normal.multiply(delta - distanceToLine));
        }
    }

    /**
     * @param point
     * @param numbersToSkip which numbers of contour part to skip
     * @return
     */
    public Point awayFromSegmentationPoint(Point point, int... numbersToSkip) {
        boolean lResult = false;
        double lAbsoluteDistance = shape.getDelta() / 2;
        double lDistance;
        double lDistance2;
        double minDistance = 10000;
        Point firstMinPoint = null;
        Point secondMinPoint = null;
        mainCycle:
        for (int i = 0; i < shape.getListOfPoints().size() - 1; i++) {
            for (int numberToSkip : numbersToSkip) {
                if (i == numberToSkip) {
                    continue mainCycle;
                }
            }
            Point firstPoint = shape.getListOfPoints().get(i);
            Point secondPoint = shape.getListOfPoints().get(i + 1);
            Point.StraightLine line = new Point.StraightLine(firstPoint, secondPoint);
            double distance = line.distanceToPoint(point);
            if (distance < minDistance) {
                minDistance = distance;
                firstMinPoint = firstPoint;
                secondMinPoint = secondPoint;
            }
        }

        if (firstMinPoint != null) {
            if (minDistance < shape.getDelta()) {
                point = awaySeg(firstMinPoint, secondMinPoint, point, point);
            }
        }

        return point;
    }

    private static class Pair<T> {
        T firstValue;
        T secondValue;

        private Pair(T firstValue, T secondValue) {
            this.firstValue = firstValue;
            this.secondValue = secondValue;
        }

        private T getFirstValue() {
            return firstValue;
        }

        private void setFirstValue(T firstValue) {
            this.firstValue = firstValue;
        }

        private T getSecondValue() {
            return secondValue;
        }

        private void setSecondValue(T secondValue) {
            this.secondValue = secondValue;
        }
    }
}
