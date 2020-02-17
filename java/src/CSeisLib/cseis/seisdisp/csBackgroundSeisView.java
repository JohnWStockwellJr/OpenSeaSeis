/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seisdisp;

import cseis.cmapgen.csCustomColorMapModel;
import cseis.general.csStandard;
import cseis.seis.csISeismicTraceBuffer;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.util.Hashtable;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * An extension of csSeisView that plots two data sets on top of each other
 * For example a seismic image with a velocity field in the background
 */
@SuppressWarnings("serial")
public class csBackgroundSeisView extends csSeisView {
  csSeisView myBackgroundSeisView;
  int myBackgroundCompositeRule;
  float myBackgroundAlpha;
  boolean myIsBackground;
  
  static final int[] COMPOSITE_RULES = {
      AlphaComposite.DST,
      AlphaComposite.DST_ATOP,
      AlphaComposite.DST_IN,
      AlphaComposite.DST_OUT,
      AlphaComposite.DST_OVER,
      AlphaComposite.SRC,
      AlphaComposite.SRC_ATOP,
      AlphaComposite.SRC_OVER,
      AlphaComposite.SRC_IN,
      AlphaComposite.SRC_OUT };
  static final String[] COMPOSITE_LABELS = {
      "AlphaComposite.DST",
      "AlphaComposite.DST_ATOP",
      "AlphaComposite.DST_IN",
      "AlphaComposite.DST_OUT",
      "AlphaComposite.DST_OVER",
      "AlphaComposite.SRC",
      "AlphaComposite.SRC_ATOP",
      "AlphaComposite.SRC_OVER",
      "AlphaComposite.SRC_IN",
      "AlphaComposite.SRC_OUT" };
  
//    csSeisView seisview = new csSeisView( this, dispSettings, myCustomColorMapModel, myColorBitType );
  public csBackgroundSeisView( JFrame parentFrame, csSeisDispSettings dispSettings, csCustomColorMapModel cmapModel, int colorMapType ) {
    super( parentFrame, dispSettings, cmapModel, colorMapType );
    myIsBackground = false;
    myBackgroundSeisView = null;
  }
  private void setParamsFromMainSeisViewStep1() {
    myBackgroundSeisView.myIsFullRepaint = (
      ( myBackgroundSeisView.myMarginLeftRight != myMarginLeftRight ) ||
      ( myBackgroundSeisView.myMarginTopBottom != myMarginTopBottom ) ||
      ( myBackgroundSeisView.mySettings.zoomHorz != mySettings.zoomHorz ) ||
      ( myBackgroundSeisView.mySettings.zoomVert != mySettings.zoomVert ) ||
      ( myBackgroundSeisView.myViewPositionHorz != myViewPositionHorz ) ||
      ( myBackgroundSeisView.myViewPositionVert != myViewPositionVert ) );

    myBackgroundSeisView.myMarginLeftRight = myMarginLeftRight;
    myBackgroundSeisView.myMarginTopBottom = myMarginTopBottom;
    myBackgroundSeisView.mySettings.zoomHorz = mySettings.zoomHorz;
    myBackgroundSeisView.mySettings.zoomVert = mySettings.zoomVert;
//    myBackgroundSeisView.setV
  }
//  private void setParamsFromMainSeisViewStep2() {
//    myBackgroundSeisView.myPrevViewPositionHorz = myPrevViewPositionHorz;
//    myBackgroundSeisView.myPrevViewPositionVert = myPrevViewPositionVert;
//    myBackgroundSeisView.myViewPositionHorz = myViewPositionHorz;
//    myBackgroundSeisView.myViewPositionVert = myViewPositionVert;
//    myBackgroundSeisView.myNextViewPositionVert = myNextViewPositionVert;
//    myBackgroundSeisView.myNextViewPositionHorz = myNextViewPositionHorz;
//  }
  private void addBackground() {
    myIsBackground = true;
    myBackgroundSeisView = new csSeisView( getParentFrame(), new csIVisibleRectListener() {
      @Override
      public Rectangle getVisibleRect() {
        return csBackgroundSeisView.this.getVisibleRect();
      }
    } );
    myBackgroundSeisView.myTEMPName = "NAME_BACKGROUND";

    myBackgroundAlpha = 0.5f;
    csSeisDispSettings ds = new csSeisDispSettings();
    ds.scaleType   = csSeisDispSettings.SCALE_TYPE_RANGE;
    ds.viType      = csSeisDispSettings.VA_TYPE_2DSPLINE;
    ds.isVIDisplay = true;
    ds.showWiggle  = false;
    ds.isNegFill   = false;
    ds.isPosFill   = false;
    myBackgroundSeisView.updateDispSettings(ds);
    myBackgroundCompositeRule = AlphaComposite.SRC_OVER;
    myBackgroundSeisView.getDispDialog().addCustomPanel( new csBackgroundDispPanel( this ) );
    // Update popup menu:
    JMenuItem itemBackgroundSettings = new JMenuItem("Display settings (background)...");
    myPopupMenu.add( itemBackgroundSettings );
    itemBackgroundSettings.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        showBackgroundSettingsDialog();
      }
    });
    myBackgroundSeisView.addSeisViewListener( new csISeisViewListener() {
      @Override
      public void changedSettings( csSeisDispSettings settings ) {
        repaint();
      }
      @Override
      public void vertScrollChanged(int scrollValue) {
      }
      @Override
      public void horzScrollChanged(int scrollValue) {
      }
      @Override
      public void sizeChanged(Dimension size) {
      }
      @Override
      public void traceBufferChanged(csISeismicTraceBuffer traceBuffer) {}
      @Override
      public void mouseExited() {}
    });
    super.addSeisViewListener( new csISeisViewListener() {
      @Override
      public void changedSettings( csSeisDispSettings settings ) {
/*        if( settings.zoomHorz != myBackgroundSeisView.mySettings.zoomHorz || settings.zoomVert != myBackgroundSeisView.mySettings.zoomVert ) {
          csSeisDispSettings newSett = new csSeisDispSettings( myBackgroundSeisView.mySettings );
          newSett.zoomHorz = settings.zoomHorz;
          newSett.zoomVert = settings.zoomVert;
          myBackgroundSeisView.updateDispSettings( newSett );
//          myBackgroundSeisView.myIsFullRepaint = true;
//          myBackgroundSeisView.repaint();
        }      */
      }
      @Override
      public void vertScrollChanged(int scrollValue) {
//        System.out.println("vert scroll changed " + scrollValue);
//        repaint();
        myBackgroundSeisView.resetViewPositionVert( scrollValue );
        myBackgroundSeisView.repaint();
      }
      @Override
      public void horzScrollChanged(int scrollValue) {
 //       System.out.println("Horz scroll changed " + scrollValue);
//        repaint();
        myBackgroundSeisView.resetViewPositionHorz( scrollValue );
        myBackgroundSeisView.repaint();
      }
      @Override
      public void sizeChanged( Dimension size ) {
        myBackgroundSeisView.setSize(size);
//        myBackgroundSeisView.repaint();
      }
      @Override
      public void traceBufferChanged( csISeismicTraceBuffer traceBuffer ) {
      }
      @Override
      public void mouseExited() {}
    });
  }
  public void updateBackground( csISeismicTraceBuffer buffer, double sampleInt ) {
    if( myBackgroundSeisView == null ) addBackground();
    myBackgroundSeisView.updateTraceBuffer( buffer, sampleInt, false );
    repaint();
  }
  public void setBackgroundCompositeRule( int compositeRule ) {
    myBackgroundCompositeRule = compositeRule;
  }
  public void setBackgroundAlpha( float alpha ) {
    myBackgroundAlpha = alpha;
  }
  public static Image makeColorTransparent( BufferedImage image, final Color color ) {
    ImageFilter filter = new RGBImageFilter() {
      // Set alpha bytes to opaque
      public int markerRGB = color.getRGB() | 0xFF000000;
      @Override
      public final int filterRGB(int x, int y, int rgb) {
        if ( ( rgb | 0xFF000000 ) == markerRGB ) {
          // Mark alpha bytes as zero (=transparent)
          return 0x00FFFFFF & rgb;
        }
        else {
          return rgb;
        }
      }
    }; 
    ImageProducer ip = new FilteredImageSource( image.getSource(), filter );
    return Toolkit.getDefaultToolkit().createImage(ip);
  }
  public void showBackgroundSettingsDialog() {
    myBackgroundSeisView.showDispSettingsDialog();
  }

  @Override
  public void paintComponent( Graphics g ) {
    Graphics2D g2 = (Graphics2D)g;
    // (0) If background data does not exist, just draw normal seismic view
//    System.out.println("(1) All views repaint");
    if( !myIsBackground ) {
        //      System.out.println("(1x) Main view repaint only");
      super.paintComponent( g );
      return;
    }
//    Rectangle rectVisible = getVisibleRect();
//    if( rectVisible.width == 0 && rectVisible.height == 0 ) return;
//    if( !myIsScrolling ) {
      setParamsFromMainSeisViewStep1();
//      System.out.println("(2) Background view repaint " + myBackgroundSeisView.myImageWidth + " " + myBackgroundSeisView.myImageHeight + " " + getVisibleRect() );
//      boolean isRepainted = myBackgroundSeisView.repaintStep1( false );
//      setParamsFromMainSeisViewStep2();

//      if( myBackgroundSeisView.myImageWidth != rectBackground.width || myBackgroundSeisView.myImageHeight != rectBackground.height ) {
//        System.out.println("(2) Background reset");
  //    if( rectBackground.width > 0 && rectBackground.height > 0 ) {
//        // (1) Paint velocity background to bitmap
//        if( !myBackgroundSeisView.resetScreenBuffer() ) return;
//        myBackgroundSeisView.myMarginLeftRight = myMarginLeftRight;
//        myBackgroundSeisView.myIsFullRepaint = true;
//        System.out.println("(3) Background full repaint");
//      }
//    if( myIsFullRepaint || myBackgroundSeisView.myIsFullRepaint ) {
//      System.out.println("(4) Background repaintStep1");
//      myBackgroundSeisView.repaintStep1( false );
      // (2) ...if seismic foreground does not include VI display:
      //     a) Paint seismic foreground (wiggle + fill) to bitmap
      //     b) Make white background of bitmap transparent
      //     c) Paint velocity background bitmap to Graphics object
      //     d) Paint seismic foreground bitmap on top
//    }
    if( !mySettings.isVIDisplay ) {
      super.paintComponent( g2 );
      Image imageSeismicWhiteTransparent = csBackgroundSeisView.makeColorTransparent( super.myBitmap, Color.white );
//      g2.drawImage(myBackgroundSeisView.myBitmap, 0, 0, this );
      g2.drawImage(myBackgroundSeisView.myVolatileBitmap, 0, 0, this );
      g2.setComposite(AlphaComposite.getInstance(myBackgroundCompositeRule, myBackgroundAlpha ) );
      g2.drawImage( imageSeismicWhiteTransparent, 0, 0, null );
    }
    // (3) ...if seismic foreground includes VI display:
    //     a) Paint velocity background to Graphics object
    //     b) Paint seismic foreground bitmap on top
    else {
        //        if( myBackgroundSeisView.myVolatileBitmap == null ) {
        //  System.out.println("(3) Drawing bitmaps " + myBackgroundSeisView.myBitmap.getWidth() + " " + myBackgroundSeisView.myBitmap.getHeight() + " NOVOLATILE ");
        // }
        // else {
        //    System.out.println("(3) Drawing bitmaps " + myBackgroundSeisView.myBitmap.getWidth() + " " + myBackgroundSeisView.myBitmap.getHeight() + " " +
        //                        myBackgroundSeisView.myVolatileBitmap.getWidth(null) + " " + 
        //                        myBackgroundSeisView.myVolatileBitmap.getHeight(null) );
        // }
      g2.drawImage(myBackgroundSeisView.myVolatileBitmap, 0, 0, this );
      g2.setComposite(AlphaComposite.getInstance(myBackgroundCompositeRule, myBackgroundAlpha ) );
      super.paintComponent( g2 );
    }
  }

  /*
  @Override
  public void zoom( float zoomVert, float zoomHorz ) {
    super.zoom(zoomVert, zoomHorz);
    if( myIsBackground ) myBackgroundSeisView.zoom(zoomVert, zoomHorz);
  }
  @Override
  public void zoom( float zoomVert, float zoomHorz, float traceIndexCentre, float sampleIndexCentre ) {
    super.zoom(zoomVert, zoomHorz, traceIndexCentre, sampleIndexCentre);
    if( myIsBackground ) myBackgroundSeisView.zoom(zoomVert, zoomHorz, traceIndexCentre, sampleIndexCentre);
  }
  @Override
  protected synchronized void resetViewPositionVert( int scrollValue ) {
    super.resetViewPositionVert( scrollValue );
    if( myIsBackground ) {
      myBackgroundSeisView.myIsScrolling = false; // Fudge
      myBackgroundSeisView.resetViewPositionVert( scrollValue );
    }
  }
  @Override
  protected synchronized void resetViewPositionHorz( int scrollValue ) {
    super.resetViewPositionHorz( scrollValue );
    if( myIsBackground ) {
      myBackgroundSeisView.myIsScrolling = false; // Fudge
      myBackgroundSeisView.resetViewPositionHorz( scrollValue );
    }
  }
  */
  //==============================================================================
  //
  public class csBackgroundDispPanel extends JPanel {
//    private final JComboBox<String> myCombo;
    private final JSlider mySliderAlpha;
    private final JTextField myTextAlpha;

    csBackgroundDispPanel( csBackgroundSeisView seisView ) {
      super( new BorderLayout() );
      setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Transparency"),
        csStandard.INNER_EMPTY_BORDER ) );

//      myCombo = new JComboBox<>( COMPOSITE_LABELS );
      myTextAlpha  = new JTextField("70");
      mySliderAlpha = new JSlider( 0, 100, 70 );
      mySliderAlpha.setMinorTickSpacing( 25 );
//    mySliderAlpha.setFont( new Font("SansSerif", Font.PLAIN, 8) );
      Hashtable<Integer,JLabel> labelTable = new Hashtable<Integer,JLabel>();
      labelTable.put( 0, new JLabel("0%") );
      labelTable.put( 50, new JLabel("50%") );
      labelTable.put( 100, new JLabel("100%") );
      mySliderAlpha.setLabelTable( labelTable );
      mySliderAlpha.setPaintLabels( true );
      mySliderAlpha.setPaintTicks( true );
//      mySliderAlpha.setPreferredSize( new Dimension(200,40) );
      myTextAlpha.setPreferredSize( new Dimension(40,myTextAlpha.getPreferredSize().height) );
      JPanel panelText = new JPanel( new GridBagLayout() );
      panelText.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, 0, 1, 1, 1.0, 0.5, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      panelText.add( myTextAlpha, new GridBagConstraints(
        0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 5, 0, 0, 5 ), 0, 0 ) );
      panelText.add( Box.createVerticalGlue(), new GridBagConstraints(
        0, 2, 1, 1, 1.0, 0.5, GridBagConstraints.WEST,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      add( panelText, BorderLayout.WEST );
      add( mySliderAlpha, BorderLayout.CENTER );
//      add( myCombo, BorderLayout.EAST );

      myTextAlpha.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          try {
            int valuePercent = Integer.parseInt( myTextAlpha.getText() );
            if( valuePercent > 100 ) valuePercent = 100;
            else if( valuePercent < 0 ) valuePercent = 0;
            myTextAlpha.setText("" + valuePercent);
            if( !mySliderAlpha.getValueIsAdjusting() ) {
              mySliderAlpha.setValue( valuePercent );
              csBackgroundSeisView.this.setBackgroundAlpha( (float)valuePercent/100.0f );
              csBackgroundSeisView.this.repaint();
            }
          }
          catch( NumberFormatException ne ) {
          }
        }
      });
      mySliderAlpha.addChangeListener( new ChangeListener() {
        @Override
        public void stateChanged(ChangeEvent e) {
          int valuePercent = mySliderAlpha.getValue();
          myTextAlpha.setText( "" + valuePercent );
          csBackgroundSeisView.this.setBackgroundAlpha( (float)valuePercent/100.0f );
          csBackgroundSeisView.this.repaint();
        }
      });
/*      myCombo.addItemListener( new ItemListener() {
        @Override
        public void itemStateChanged(ItemEvent e) {
          if( e.getStateChange() == ItemEvent.SELECTED ) {
            int index = myCombo.getSelectedIndex();
            if( index >= 0 && index < COMPOSITE_RULES.length ) {
              csBackgroundSeisView.this.setBackgroundCompositeRule( COMPOSITE_RULES[index] );
              csBackgroundSeisView.this.repaint();
            }
          }
        }
      }); */
    }
  }

}
