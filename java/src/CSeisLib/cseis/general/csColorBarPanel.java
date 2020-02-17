/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import javax.swing.*;
import java.text.DecimalFormat;
import java.awt.*;

/**
 * Color bar.<br>
 * Plots color bar as a JPanel.
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csColorBarPanel extends JPanel {
  public static final int ORIENT_HORIZONTAL   = 1;
  public static final int ORIENT_VERTICAL     = 2;
  public static final int ANNOTATION_NONE     = -1;
  public static final int ANNOTATION_SIMPLE   = 3;
  public static final int ANNOTATION_ADVANCED = 4;

  public static final int NUM_COLOR_RESOLUTION = 1001;
  
  private csCustomColorMap myMap;
  private int myOrientation;
  private int myAnnotation;
  private int myFontSize = 12;
  private double myValueStep;
  private boolean myIsAutoSize;
  private boolean myShowKneePoints;

  public csColorBarPanel( csColorMap map ) {
    this( map, csColorBarPanel.ORIENT_VERTICAL );
  }
  public csColorBarPanel( csColorMap map, int orientation ) {
    this( map, csColorBarPanel.ORIENT_VERTICAL, ANNOTATION_SIMPLE );
  }
  public csColorBarPanel( csColorMap map, int orientation, int annotation ) {
    super( new BorderLayout() );
    myMap = new csCustomColorMap( map );
    myMap.setNumColorResolution( NUM_COLOR_RESOLUTION );
    myOrientation = orientation;
    myAnnotation  = annotation;
    myIsAutoSize = false;
    myShowKneePoints = false;
    setPreferredSize( new Dimension(50,0) );
    setMinimumSize(new Dimension(0,0));
  }
  public void setAnnotation( int annotation ) {
    myAnnotation = annotation;
  }
  public void setShowKneePoints( boolean doShow ) {
    myShowKneePoints = doShow;
  }
  public void setFontSize( int fontSize ) {
    myFontSize = fontSize;
  }
  public void setMinMax( double minValue, double maxValue ) {
    myMap.setMinMax( minValue, maxValue );
  }
  public csCustomColorMap getColorMap() {
    return myMap;
  }
  public void updateColors( Color[] colors, double[] weights ) {
    myMap.setColors( colors, weights, myMap.getNumColorResolution() );
    revalidate();
    repaint();
  }
  public void setColorMap( csColorMap map ) {
    if( map.getMinValue() == myMap.getMinValue() && map.getMaxValue() == myMap.getMaxValue() ) {
      Color[] colors1 = myMap.getColorKneePoints();
      Color[] colors2 = map.getColorKneePoints();
      if( colors1.length == colors2.length ) {
        boolean isEqual = true;
        for( int i = 0; i < colors1.length; i++ ) {
          if( colors1[i] != colors2[i] || myMap.getWeightKneePoint(i) != map.getWeightKneePoint(i) ) {
            isEqual = false;
            break;
          }
        }
        if( isEqual ) return;
      }
    }
    myMap = new csCustomColorMap( map );
    myMap.setNumColorResolution( NUM_COLOR_RESOLUTION );
    revalidate();
    repaint();
  }
  /**
   * Set auto-size to resize this component automatically when it is too small to fit into its container
   * @param doAutoSize 
   */
  public void setAutoSize( boolean doAutoSize ) {
    myIsAutoSize = doAutoSize;
  }
  public void paintComponent( Graphics g1 ) {
    super.paintComponent(g1);
    double minValue = myMap.getMinValue();
    double maxValue = myMap.getMaxValue();

    Graphics2D g = (Graphics2D)g1;
    g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
    g.setStroke( new BasicStroke(1.0f) );
    g.setFont(new Font("SansSerif", Font.PLAIN, myFontSize));
    g.setColor(Color.white);

    DecimalFormat formatter = setAnnFormatter( maxValue-minValue );
    // Re-evaluate component width. Enlarge if needed:
    if( myIsAutoSize && myAnnotation == ANNOTATION_ADVANCED ) {
      String textMin = formatter.format( minValue );
      String textMax = formatter.format( maxValue );
      int fontWidth  = Math.max( g1.getFontMetrics().stringWidth( textMin ), g1.getFontMetrics().stringWidth( textMax ) );
      int idealWidth = 50 + fontWidth;
      if( myOrientation == ORIENT_VERTICAL && idealWidth > getWidth() ) {
        setPreferredSize( new Dimension( idealWidth, getPreferredSize().height ) );
      }
    }
    
    int dist = 10;
    int height = this.getHeight();
    int width = this.getWidth();
    g.fillRect(0, 0, width, height);

    if( myAnnotation == ANNOTATION_ADVANCED ) {
      // Draw black box
      g.setColor( Color.black );
      g.drawRect( dist-1, dist-1, width/3+2, height-(2*dist)+2);
  
      int minY = dist;
      int maxY = height-dist;
      int diff = maxY - minY;
      if( diff > 0 ) {
        for( int yPos = minY; yPos <= maxY; yPos++ ) {
          float s = (float)(yPos - minY) / (float)diff;
          float value = (float)(maxValue + s*( minValue - maxValue ));
          g.setColor( myMap.getColor( value ) );
          g.drawLine( dist, yPos, dist+width/3, yPos);
        }
        if( myShowKneePoints ) {
          drawKneePoints( g, minValue, maxValue, diff, dist, minY );
        }
/*        if( myShowKneePoints ) {
          Stroke savedStroke = g.getStroke();
          g.setStroke( new BasicStroke(3.0f) );
          for( int ip = 1; ip < myMap.getNumKneePoints()-1; ip++ ) {
            float pos = (float)myMap.getWeightKneePoint(ip);
            // Convert pos into yPos:
            float s = (float)( (pos - maxValue) / (minValue - maxValue) );
            int yPos = (int)Math.round( s * diff ) + minY + 1;  // +1 to compensate for stroke size of 3
            g.setColor( myMap.getColorKneePoint(ip) );
            g.drawLine( 0, yPos, dist-1, yPos);
            g.setColor( Color.black );
            g.drawLine( dist, yPos, 2*dist, yPos);
          }
          g.setStroke( savedStroke );
        } */
      }
      int labelInc = 120;
      int maxSteps = (int)( (float)(maxY-minY)/(float)(labelInc) );
      int nn = (int)((maxValue-minValue)/((float)maxSteps*myValueStep));
      if( nn <= 0 ) nn = 1;
      myValueStep = nn*myValueStep;
      double yStep = (double)(maxY-minY)/(maxValue-minValue)*myValueStep-0.3;
      int nSteps = (int)((maxValue-minValue)/myValueStep + 0.5) + 1;
      if( nSteps < 0 ) nSteps = 0;

      g.setColor(Color.black);
      int fontHeight = g.getFontMetrics().getHeight();
      for( int istep = 0; istep < nSteps; istep++ ) {
        int yPos = (int)(-istep*yStep) + maxY+1;
        double value = minValue + istep*myValueStep;
        if( istep == nSteps-1 ) {
          yPos = minY-1;
          value = maxValue;
        }
        String text = formatter.format( value );
        g.drawLine( dist+width/3+2, yPos, (dist+width/3+2+5), yPos);
        if( istep == 0 ) yPos -= fontHeight/3;
        g.drawString( text, 2*dist+width/3, yPos+(int)(fontHeight/2.5) );
      }
    }
    else { // if( myAnnotation == ANNOTATION_SIMPLE ) {
      // Draw black box
      g.setColor( Color.black );
      g.drawRect( dist-1, 2*dist-1, width-2*(dist-1), height-4*dist+2);
  
      int minY = 2*dist;
      int maxY = height-2*dist;
      int diff = maxY - minY + 1;
      if( diff > 0 ) {
        for( int yPos = minY; yPos <= maxY; yPos++ ) {
          float s = (float)(yPos - minY) / (float)diff;
          float value = (float)(maxValue + s*( minValue - maxValue ));
          g.setColor( myMap.getColor( value ) );
          g.drawLine( dist, yPos, width-dist, yPos);
        }
        if( myShowKneePoints ) {
          drawKneePoints( g, minValue, maxValue, diff, dist, minY );
        }
      }
      
      g.setColor(Color.black);
      int fontHeight = g.getFontMetrics().getHeight();

      String textMin = formatter.format( minValue );
      String textMax = formatter.format( maxValue );

      int yPos = dist;
      g.drawString( textMax, dist, yPos+(int)(fontHeight/2.5) );
      yPos = height - dist;
      g.drawString( textMin, dist, yPos+(int)(fontHeight/2.5) );
    }
  }
  private void drawKneePoints( Graphics2D g, double minValue, double maxValue, double diff, int dist, int minY ) {
    Stroke savedStroke = g.getStroke();
    g.setStroke( new BasicStroke(3.0f) );
    for( int ip = 1; ip < myMap.getNumKneePoints()-1; ip++ ) {
      float pos = (float)myMap.getWeightKneePoint(ip);
      // Convert pos into yPos:
      float s = (float)( (pos - maxValue) / (minValue - maxValue) );
      int yPos = (int)Math.round( s * diff ) + minY + 1;  // +1 to compensate for stroke size of 3
      g.setColor( myMap.getColorKneePoint(ip) );
      g.drawLine( 0, yPos, dist-1, yPos);
      g.setColor( Color.black );
      g.drawLine( dist, yPos, 2*dist, yPos);
    }
    g.setStroke( savedStroke );
  }
  private DecimalFormat setAnnFormatter( double valueStep ) {
    String format = "";
    if( myMap.getSmoothMode() == csColorMap.SMOOTH_MODE_DISCRETE ) {
      valueStep = (double)myMap.getDiscreteStep();
      format = "0";
    }
    else if( valueStep > 100 ) {
      valueStep = 10;
      format = "0";
    }
    else if( valueStep > 10 ) {
      valueStep = 0.5;
      format = "0.0";
    }
    else if( valueStep > 1 ) {
      valueStep = 0.05;
      format = "0.00";
    }
    else if( valueStep > 0.1 ) {
      valueStep = 0.005;
      format = "0.000";
    }
    else if( valueStep > 0.01 ) {
      valueStep = 0.0005;
      format = "0.0000";
    }
    else if( valueStep > 0.001 ) {
      format = "0.00000";
      valueStep = 0.00005;
    }
    else {
      format = "";
      valueStep /= 10.0;
    }
    DecimalFormat formatter = new DecimalFormat(format);
    myValueStep = valueStep;
    return formatter;
  }
}


