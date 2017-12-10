/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
//import static javax.swing.Spring.scale;
//import static jogamp.nativewindow.SurfaceScaleUtils.scale;
import util.TFChangeListener;
import util.VectorMath;

import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    
    private boolean shading = false;  // color of the light
    private TFColor I_a = new TFColor();
    private TFColor I_d = new TFColor();
    private int method = 1;
    private double k_amb = 0.1;//0.1;
    private double k_dif = 0.7;////0.7;
    private double k_spec = 0.2;
    private double alpha1 = 10;
    private boolean isTrilinearInterpolation = false;
 
    public void setShading()
    {
     if(shading == true){shading = false; this.changed();}
     else{
     shading = true; 
     this.changed();
     }
    }
    public void setTrilinear()
    {
     if(isTrilinearInterpolation == true){isTrilinearInterpolation = false; this.changed();}
     else{
     isTrilinearInterpolation = true; 
     this.changed();
     }
    }

    // Blinn-Phong shading model

    
    public TFColor getColor(double[] L, double[] N, double[] V, TFColor voxelColor) {
        // compute reflection vector
        // R = 2(N.L)N - L
        double[] R = VectorMath.subtract((VectorMath.scale(N, 2*VectorMath.dotproduct(N, L))), L);

        // compute L.N and (V.R)^{alpha}
        double L_dot_N = Math.max(0, VectorMath.dotproduct(L, N));
        double V_dot_R = VectorMath.dotproduct(V, R);
        double V_dot_R_power_alpha = (V_dot_R <= 0 || L_dot_N <= 0) ? 0 : Math.pow(V_dot_R, alpha1);
        I_d = voxelColor;
        I_a = new TFColor(1,1,1,1);
        // compute color
        TFColor color = new TFColor(
                k_amb * I_a.r + k_dif * I_d.r * L_dot_N + k_spec * I_d.r * V_dot_R_power_alpha,
                k_amb * I_a.g + k_dif * I_d.g * L_dot_N + k_spec * I_d.g * V_dot_R_power_alpha,
                k_amb * I_a.b + k_dif * I_d.b * L_dot_N + k_spec * I_d.b * V_dot_R_power_alpha,
                I_d.a);
        return color;
    }
    private TFColor getShadedColor(double[] coord, TFColor voxelColor, double[] viewVec) {
        if (coord[0] < 0 || coord[0] >= volume.getDimX()
                || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            voxelColor.set(0, 0, 0, voxelColor.a);
            return voxelColor;
        }
        I_d = voxelColor;
        VoxelGradient vg = gradients.getGradient(
                (int) Math.floor(coord[0]),
                (int) Math.floor(coord[1]),
                (int) Math.floor(coord[2]));
        return getColor(viewVec, vg.getNormal(), viewVec, I_d);
        
    }
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    private double getGradient(double[] coord) {
        if (coord[0] < 0 || coord[0] >= volume.getDimX() - 1
                || coord[1] < 0 || coord[1] >= volume.getDimY() - 1
                || coord[2] < 0 || coord[2] >= volume.getDimZ() - 1) {
            return 0;
        }

        if (isTrilinearInterpolation) {
            int x1 = (int) Math.floor(coord[0]);
            int y1 = (int) Math.floor(coord[1]);
            int z1 = (int) Math.floor(coord[2]);
            int x2 = x1 + 1;
            int y2 = y1 + 1;
            int z2 = z1 + 1;

            double alpha = coord[0] - x1;
            double beta = coord[1] - y1;
            double gamma = coord[2] - z1;

            return (1 - alpha) * (1 - beta) * (1 - gamma) * gradients.getGradient(x1, y1, z1).mag +
                    alpha * (1 - beta) * (1 - gamma) * gradients.getGradient(x2, y1, z1).mag +
                    (1 - alpha) * beta * (1 - gamma) * gradients.getGradient(x1, y2, z1).mag +
                    alpha * beta * (1 - gamma) * gradients.getGradient(x2, y2, z1).mag +
                    (1 - alpha) * (1 - beta) * gamma * gradients.getGradient(x1, y1, z2).mag +
                    alpha * (1 - beta) * gamma * gradients.getGradient(x2, y1, z2).mag +
                    (1 - alpha) * beta * gamma * gradients.getGradient(x1, y2, z2).mag +
                    alpha * beta * gamma * gradients.getGradient(x2, y2, z2).mag;
        } else {
            int x = (int) Math.floor(coord[0]);
            int y = (int) Math.floor(coord[1]);
            int z = (int) Math.floor(coord[2]);
            return gradients.getGradient(x, y, z).mag;
        }
    }
       
    double getVoxel(double[] coord) {

        // x y z in data
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        
        // x y z Floored and Ceiled 
        int x0 = (int) Math.floor(x);
        int x1 = (int) Math.min(Math.ceil(x), volume.getDimX() - 1);
        int y0 = (int) Math.floor(y);
        int y1 = (int) Math.min(Math.ceil(y), volume.getDimY() - 1);
        int z0 = (int) Math.floor(z);
        int z1 = (int) Math.min(Math.ceil(z), volume.getDimZ() - 1);
        
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ() || x0 >= volume.getDimX() || x0 < 0 || y0 >= volume.getDimY() || y0 < 0 || z0 >= volume.getDimZ() || z0<0 || x1<0 || x1 > volume.getDimX() || y1 < 0 || y1 > volume.getDimY() || z1<0 || z1> volume.getDimZ()) {
            return 0;
        }
        
        if(isTrilinearInterpolation){
            // If trilinear implementation is chosen
        
            double sx_000 = volume.getVoxel(x0, y0, z0);
            double sx_001 = volume.getVoxel(x0, y0, z1);
            double sx_010 = volume.getVoxel(x0, y1, z0);
            double sx_100 = volume.getVoxel(x1, y0, z0);
            double sx_101 = volume.getVoxel(x1, y0, z1);
            double sx_011 = volume.getVoxel(x0, y1, z1);
            double sx_110 = volume.getVoxel(x1, y1, z0);
            double sx_111 = volume.getVoxel(x1, y1, z1);

            // interpolate along x dimension
            double sx_00 = linearInterpolate(x, x0, x1, sx_000, sx_100);
            double sx_01 = linearInterpolate(x, x0, x1, sx_001, sx_101);
            double sx_10 = linearInterpolate(x, x0, x1, sx_010, sx_110);
            double sx_11 = linearInterpolate(x, x0, x1, sx_011, sx_111);

            // interpolate along y dimension
            double sx_0 = linearInterpolate(y, y0, y1, sx_00, sx_10);
            double sx_1 = linearInterpolate(y, y0, y1, sx_01, sx_11);

            // interpolate along z dimension
            double sx = linearInterpolate(z, z0, z1, sx_0, sx_1);
            return sx;
        }
        else{
            return volume.getVoxel(x0, y0, z0);
        }
    }
    
    public static double linearInterpolate(double x, double x0, double x1, double v0, double v1){
           double alpha = (x-x0)/(x1-x0);
           return (1 - alpha)*v0 + alpha*v1;
    }  

    void slicer(double[] viewMatrix) {
        int resolution = 1;
        if(interactiveMode) {resolution = 3;}
        
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j=j+resolution) {
            for (int i = 0; i < image.getWidth(); i=i+resolution) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                double val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                for (int m = 0; m <resolution; m++) {
                    for (int n= 0; n <resolution; n++) {
                        if (n + i < image.getHeight() 
                            && n + i >= 0
                            && m + j < image.getWidth() 
                            && m + j >= 0) {
                            image.setRGB(n + i, m + j,pixelColor);
                        }
                    }
                }
            }
        }

    }
    


    void maximumIntensityProjection(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        
        int resolution = 1;
        if (interactiveMode) {
            resolution = 1  ;
        }
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] viewVec = new double[3];     
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        for (int j = 0; j < image.getHeight(); j+= resolution) {
            for (int i = 0; i < image.getWidth(); i+= resolution) {
                
                int diag = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()); 
                double maxIntensity = 0;
                
                for (int k = 0; k < diag - 1; k++) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) +(k - imageCenter) * viewVec[0] + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + (k - imageCenter) * viewVec[1] + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) +  (k - imageCenter)*viewVec[2] + volumeCenter[2];


                    double val = getVoxel(pixelCoord);
                
                
                    if (val > maxIntensity) {
                        maxIntensity = val;
                    }
                    if (val / max > 0.95) {
                        break;
                    }
                }
               
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxIntensity/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxIntensity > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                
                for (int m = 0; m <resolution; m++) {
                    for (int n= 0; n <resolution; n++) {
                        if (n + i < image.getHeight() 
                            && n + i >= 0
                            && m + j < image.getWidth() 
                            && m + j >= 0) {
                            image.setRGB(n + i, m + j,pixelColor);
                        }
                    }
                }
            }
        }

    }

    private void compositing(double[] viewMatrix) {
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] viewVec = new double[3];
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.normalize(viewVec);
        
        int resolution = 1;
        if(interactiveMode) {
            resolution = 3;
        }
        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        TFColor voxelColor, voxelColor1 = new TFColor();
        

        // ray starting coordinates
        double[] rayCoord = new double[3];

        int diag = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ());
 
        //float diagLen = (float)volume.getDiagonalLength();
        for (int j = 0; j < image.getHeight(); j=j+resolution) {
            for (int i = 0; i < image.getWidth(); i=i+resolution) {
                // Create a ray


                // Check intersection points with the ray
               

                voxelColor1.set(0, 0, 0, 1);
                //double[] v = VectorMath.subtract(farP, nearP);
                 
                
                for (int step=0;step < diag - 1;step = step + 1) {
                rayCoord[0] = uVec[0] * (i - imageCenter) 
                        + vVec[0] * (j - imageCenter) 
                         +  viewVec[0]*(step - imageCenter ) + volumeCenter[0];
                rayCoord[1] = uVec[1] * (i - imageCenter)
                        + vVec[1] * (j - imageCenter) 
                         + viewVec[1]*(step - imageCenter) + volumeCenter[1];
                rayCoord[2] = uVec[2] * (i - imageCenter) 
                        + vVec[2] * (j - imageCenter)
                         + viewVec[2]*(step - imageCenter)+volumeCenter[2];

                    double val = getVoxel(rayCoord);

                    // apply the transfer function to obtain a color
                    voxelColor = tFunc.getColor(val);

                    // phong shading
                    if (shading) {
                        voxelColor = getShadedColor(rayCoord, voxelColor, viewVec);
                    }
                    
                    // compositing
                    voxelColor1.r = voxelColor.a * voxelColor.r + (1 - voxelColor.a) * voxelColor1.r;
                    voxelColor1.g = voxelColor.a * voxelColor.g + (1 - voxelColor.a) * voxelColor1.g;
                    voxelColor1.b = voxelColor.a * voxelColor.b + (1 - voxelColor.a) * voxelColor1.b;
                }
 
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor1.a <= 1.0 ? (int) Math.floor(voxelColor1.a * 255) : 255;
                int c_red = voxelColor1.r <= 1.0 ? (int) Math.floor(voxelColor1.r * 255) : 255;
                int c_green = voxelColor1.g <= 1.0 ? (int) Math.floor(voxelColor1.g * 255) : 255;
                int c_blue = voxelColor1.b <= 1.0 ? (int) Math.floor(voxelColor1.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                for (int m = 0; m <resolution; m++) {
                    for (int n= 0; n <resolution; n++) {
                        if (n + i < image.getHeight() 
                            && n + i >= 0
                            && m + j < image.getWidth() 
                            && m + j >= 0) {
                            image.setRGB(n + i, m + j,pixelColor);
                        }
                    }
                }
            }
        }
    }

     void transferTwoD (double[] viewMatrix){
          // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        int resolution = 1;
        if (interactiveMode) {
            resolution = 3  ;
        }
        
                // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];  
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
    //    TFColor voxelColor = new TFColor();

        int diag = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ());

        // For each pixel of the image
        for (int j = 0; j < image.getHeight(); j+=resolution) {
            for (int i = 0; i < image.getWidth(); i+=resolution) {
                // First, set a color variable in which we can put all the colors together
                TFColor voxelColor = new TFColor(0,0,0,1);
                
                for(int k = 0; k < diag - 1; k++) {
                    // Calculate the coordinates of the pixel in the data
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (k - imageCenter) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (k - imageCenter) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (k - imageCenter) + volumeCenter[2];
                  
                    int val = (int) getVoxel(pixelCoord);
                   
                    if (val == 0) {
                        continue;
                    }
                    float baseIntensity = tfEditor2D.triangleWidget.baseIntensity;
                    TFColor color = tfEditor2D.triangleWidget.color;
                    double radius = tfEditor2D.triangleWidget.radius;
                    double maxGradient = tfEditor2D.triangleWidget.maxGradient;
                    double minGradient = tfEditor2D.triangleWidget.minGradient;
                    double alpha;
                    
                    VoxelGradient gradient = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);
                    if (gradient.mag > maxGradient || gradient.mag < minGradient)
                       {alpha=0;}
                    else if (val == baseIntensity && gradient.mag == 0){
                        alpha = 1;
                    } else if (gradient.mag > 0 && ((val - (radius * gradient.mag)) <= baseIntensity && baseIntensity <= (val + (radius * gradient.mag)))) {
                            alpha = color.a * (1 - (1 / radius) * Math.abs(((baseIntensity - val) / gradient.mag)));
                        } else {
                        alpha = 0;
                    }
                    
                    
                    if (shading) {
                        color = getShadedColor(pixelCoord, color,viewVec);
                    } else {
                        color= color;
                    }
                    voxelColor.r = (color.r * alpha ) + (voxelColor.r * (1 - alpha));
                    voxelColor.g = (color.g * alpha ) + (voxelColor.g * (1 - alpha));
                    voxelColor.b = (color.b * alpha ) + (voxelColor.b * (1 - alpha));    
                }
                
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                
               for (int m = 0; m <resolution; m++) {
                    for (int n= 0; n <resolution; n++) {
                        if (n + i < image.getHeight() 
                            && n + i >= 0
                            && m + j < image.getWidth() 
                            && m + j >= 0) {
                            image.setRGB(n + i, m + j,pixelColor);
                        }
                    }
                }
            }
        }
     }
  

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
            
            switch (method){
            case 1:
                slicer(viewMatrix);
                break;
            case 2:
                maximumIntensityProjection(viewMatrix);
                break;
            case 3:
                compositing(viewMatrix);
                break;
            case 4:
                transferTwoD(viewMatrix);
                break;
            default:
                slicer(viewMatrix);
        }    
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];
    public void setMethod(int method) {
        this.method = method;
        this.changed();
    }
    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
