/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.general;

import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.JButton;
import java.awt.Color;
import java.awt.image.BufferedImage;
import javax.swing.ImageIcon;

/**
 * Standard CSeisLib parameters.
 * @author Bjorn Olofsson
 */
public class csStandard {
  public static final Border INNER_EMPTY_BORDER = BorderFactory.createEmptyBorder( 1, 3, 3, 3 );
  public static final Border DIALOG_BORDER = BorderFactory.createEmptyBorder( 4, 7, 7, 7 );
  public static final Border INNER_PARAM_BORDER = BorderFactory.createCompoundBorder(
      BorderFactory.createCompoundBorder(BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder()),
      csStandard.DIALOG_BORDER );

  public static JButton createColorButton( Color color ) {
    ImageIcon icon = new ImageIcon( new BufferedImage( 12, 12, BufferedImage.TYPE_INT_RGB ) );
    JButton b = new JButton( icon );
    return b;
  }
  
 public static final double ABSENT_VALUE      = Double.MAX_VALUE;
 public static final float ABSENT_VALUE_FLOAT = Float.MAX_VALUE;
 public static final int ABSENT_VALUE_INT     = Integer.MAX_VALUE;
 public static final String ABSENT_VALUE_STR  = "NaN";

 
 public static final int ICON_SIZE = 16;

}


