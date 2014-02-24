package shape;

import draw.Color;
import draw.ColorManager;
import draw.Helper;
import processing.core.PApplet;

import java.text.DecimalFormat;

/**
 * User: shyyko
 * Date: 24.02.14
 * Time: 19:06
 */
public class ColorManagerDrawer {
    public static final int SQUARE_SIZE = 20;
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.##");

    private int top;
    private Helper helper;
    private PApplet drawer;
    private ColorManager colorManager;
    private int left;

    public ColorManagerDrawer(Helper helper, PApplet drawer, ColorManager colorManager, int top, int left) {
        this.helper = helper;
        this.drawer = drawer;
        this.colorManager = colorManager;
        this.top = top;
        this.left = left;
    }

    public void setColorManager(ColorManager colorManager) {
        this.colorManager = colorManager;
    }

    public void draw() {
        int top = this.top;
        Color color = colorManager.getMinimumColor();
        double value = colorManager.getMinimalValue();
        for (int i = 0; i < ColorManager.COLOR_STEPS; i++) {
            drawer.fill(color.getR(), color.getG(), color.getB());
            drawer.rect(left, top, SQUARE_SIZE, SQUARE_SIZE);
            drawer.text(DECIMAL_FORMAT.format(value), left + SQUARE_SIZE * 1.5f, top);
            color = color.plus(colorManager.getColorStep());
            top += SQUARE_SIZE;
            value += colorManager.getValueStep();
        }
        drawer.text(DECIMAL_FORMAT.format(value), left + SQUARE_SIZE * 1.5f, top);
    }
}
