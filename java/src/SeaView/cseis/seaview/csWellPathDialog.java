package cseis.seaview;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import cseis.math.csPoint3D;
import cseis.seisdisp.csWellPathAttr;
import cseis.seisdisp.csWellPathOverlayAttr;
import cseis.swing.csColorButton;
import cseis.swing.csColorChangeListener;


@SuppressWarnings("serial")
public class csWellPathDialog extends JDialog {
  public static final int MIN_WELL_FILES = 1;
  private static int BUNDLE_ID_COUNTER = 0;
  static private Color[] COLOR_LIST = { Color.red, Color.blue, Color.green, Color.orange, Color.yellow, Color.black, Color.white };

  private ArrayList<WellPathBundle> myWellPathBundleList;
  private JPanel myMainPanel;
  private JButton myButtonAdd;
  private JButton myButtonClose;
  private JButton myButtonApply;
  private JCheckBox myBoxShowShadow;
  private JSlider mySliderShadowDepth;

  private JTextField myTextHdrX;
  private JTextField myTextHdrY;
  private JTextField myTextLineSize;
  private JTextField myTextMaxDistance;
  private JSlider mySliderTransparency;

  private String myCurrentDirectory;
  private ArrayList<csIWellPathDialogListener> myListeners;
  private JFrame myParentFrame;
  
  public csWellPathDialog( JFrame parentFrame, String directory ) {
    super( parentFrame, "Well path setup" );
    myParentFrame = parentFrame;
    myWellPathBundleList = new ArrayList<WellPathBundle>(MIN_WELL_FILES);
    myMainPanel = new JPanel( new GridBagLayout() );
    myCurrentDirectory = directory;
    myListeners = new ArrayList<csIWellPathDialogListener>(1);
    for( int i = 0; i < MIN_WELL_FILES; i++ ) {
      myWellPathBundleList.add( new WellPathBundle() );
    }

    myTextHdrX = new JTextField("bin_x");
    myTextHdrY = new JTextField("bin_y");
    myTextLineSize = new JTextField("3");
    myTextMaxDistance = new JTextField("5000");
    myTextHdrX.setToolTipText("Name of trace header containing X coordinate");
    myTextHdrY.setToolTipText("Name of trace header containing Y coordinate");
    myTextLineSize.setToolTipText("Width of plotted well path (in pixels)");
    myTextMaxDistance.setToolTipText("<html>Maximum distance from seismic section to plot.<br>" + 
        "<i>Only the well path portion which is close enough to the current seismic section will be displayed</i></html>");

    myButtonAdd   = new JButton("Add");
    myButtonApply = new JButton("Apply");
    myButtonClose = new JButton("Close");
    myBoxShowShadow = new JCheckBox("Shadow");
    myButtonAdd.setToolTipText( "Add well path" );
    myButtonApply.setToolTipText( "Apply changes" );
    myButtonClose.setToolTipText( "Close dialog window" );
    myBoxShowShadow.setToolTipText( "Show shadow from well path" );

    mySliderShadowDepth= new JSlider( JSlider.HORIZONTAL, csWellPathOverlayAttr.MIN_SHADOW_DEPTH, csWellPathOverlayAttr.MAX_SHADOW_DEPTH,
        new csWellPathOverlayAttr().shadowDepth );
    mySliderShadowDepth.setToolTipText("Shadow depth");
    mySliderShadowDepth.setMajorTickSpacing( csWellPathOverlayAttr.MAX_SHADOW_DEPTH - csWellPathOverlayAttr.MIN_SHADOW_DEPTH );
    mySliderShadowDepth.setPaintTicks(true);
    mySliderShadowDepth.setToolTipText("Distance that shadow is cast away from well path");

    Hashtable<Integer,JLabel> labelsSlider1 = new Hashtable<Integer,JLabel>();
    labelsSlider1.put( csWellPathOverlayAttr.MIN_SHADOW_DEPTH, new JLabel("Short") );
    labelsSlider1.put( csWellPathOverlayAttr.MAX_SHADOW_DEPTH, new JLabel("Long") );
    mySliderShadowDepth.setLabelTable( labelsSlider1 );
    mySliderShadowDepth.setPaintLabels(true);
    mySliderShadowDepth.setSnapToTicks(false);
    
    mySliderTransparency = new JSlider( JSlider.HORIZONTAL, csWellPathOverlayAttr.MIN_TRANSPARENCY, csWellPathOverlayAttr.MAX_TRANSPARENCY,
        new csWellPathOverlayAttr().transparency );
    mySliderTransparency.setToolTipText("Transparency of well path plot");
    mySliderTransparency.setMajorTickSpacing( csWellPathOverlayAttr.MAX_TRANSPARENCY - csWellPathOverlayAttr.MIN_TRANSPARENCY );
    mySliderTransparency.setPaintTicks(true);
    mySliderTransparency.setToolTipText("Brightness/opacity of well path plot");

    Hashtable<Integer,JLabel> labelsSlider = new Hashtable<Integer,JLabel>();
    labelsSlider.put( csWellPathOverlayAttr.MIN_TRANSPARENCY, new JLabel("Dim") );
    labelsSlider.put( csWellPathOverlayAttr.MAX_TRANSPARENCY, new JLabel("Bright") );
    mySliderTransparency.setLabelTable( labelsSlider );
    mySliderTransparency.setPaintLabels(true);
    mySliderTransparency.setSnapToTicks(false);
    
    updatePanel();

    JPanel panelHdr = new JPanel( new GridBagLayout() );
    panelHdr.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createEmptyBorder(5,10,5,10),
        BorderFactory.createTitledBorder("Trace header setup")));
 
    JPanel panelPlotParam = new JPanel( new GridBagLayout() );
    panelPlotParam.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createEmptyBorder(5,10,5,10),
        BorderFactory.createTitledBorder("Plotting parameters")));
    
    JPanel panelCenter= new JPanel( new BorderLayout() );
    panelCenter.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createEmptyBorder(5,10,5,10),
        BorderFactory.createTitledBorder("Well path file setup")));

    panelHdr.add( new JLabel("X coordinate:"), new GridBagConstraints(
        0,0,1,1,0.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,5,0,10),0,0));
    panelHdr.add( myTextHdrX, new GridBagConstraints(
        1,0,1,1,1.0,0,GridBagConstraints.WEST,GridBagConstraints.HORIZONTAL,new Insets(5,0,0,5),0,0));
    panelHdr.add( new JLabel("Y coordinate:"), new GridBagConstraints(
        0,1,1,1,0.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,5,5,10),0,0));
    panelHdr.add( myTextHdrY, new GridBagConstraints(
        1,1,1,1,1.0,0,GridBagConstraints.WEST,GridBagConstraints.HORIZONTAL,new Insets(5,0,5,5),0,0));

    int yp = 0;
    panelPlotParam.add( new JLabel("Line width:"), new GridBagConstraints(
        0,0,1,1,0.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,5,0,10),0,0));
    panelPlotParam.add( myTextLineSize, new GridBagConstraints(
        1,yp++,1,1,1.0,0,GridBagConstraints.WEST,GridBagConstraints.HORIZONTAL,new Insets(5,0,0,5),0,0));
    panelPlotParam.add( new JLabel("Max distance [m]:"), new GridBagConstraints(
        0,yp,1,1,0.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,5,5,10),0,0));
    panelPlotParam.add( myTextMaxDistance, new GridBagConstraints(
        1,yp++,1,1,1.0,0,GridBagConstraints.WEST,GridBagConstraints.HORIZONTAL,new Insets(5,0,5,5),0,0));
    panelPlotParam.add( new JLabel("Transparency:"), new GridBagConstraints(
        0,yp,1,1,0.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,5,5,10),0,0));
    panelPlotParam.add( mySliderTransparency, new GridBagConstraints(
        1,yp++,1,1,1.0,0,GridBagConstraints.WEST,GridBagConstraints.HORIZONTAL,new Insets(5,0,5,5),0,0));
    panelPlotParam.add( myBoxShowShadow, new GridBagConstraints(
        0,yp,1,1,0.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,5,5,10),0,0));
    panelPlotParam.add( mySliderShadowDepth, new GridBagConstraints(
        1,yp++,1,1,1.0,0,GridBagConstraints.WEST,GridBagConstraints.HORIZONTAL,new Insets(5,0,5,5),0,0));

    JScrollPane pane = new JScrollPane( myMainPanel, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
    panelCenter.add( pane );

    JPanel panelTop = new JPanel( new BorderLayout() );
    panelTop.add( panelHdr, BorderLayout.NORTH );
    panelTop.add( panelPlotParam, BorderLayout.CENTER);
    
    JPanel panelButtons = new JPanel( new GridBagLayout() );
    int xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++,0,1,1,1.0,0,GridBagConstraints.WEST,GridBagConstraints.BOTH,new Insets(5,0,5,0),0,0));
    panelButtons.add( myButtonApply, new GridBagConstraints(
        xp++,0,1,1,0.0,0,GridBagConstraints.EAST,GridBagConstraints.NONE,new Insets(5,0,5,10),0,0));
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++,0,1,1,0.0,0,GridBagConstraints.EAST,GridBagConstraints.NONE,new Insets(5,0,5,10),0,0));
    
    getContentPane().add( panelTop, BorderLayout.NORTH );
    getContentPane().add( panelCenter, BorderLayout.CENTER );
    getContentPane().add( panelButtons, BorderLayout.SOUTH );
    pack();
    setLocationRelativeTo( parentFrame );
    
    myButtonAdd.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        addBundle();
      }
    });
    myButtonApply.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    myButtonClose.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        csWellPathDialog.this.dispose();
      }
    });
    myBoxShowShadow.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    mySliderShadowDepth.addChangeListener( new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        if( !mySliderShadowDepth.getValueIsAdjusting() ) {
          apply();
        }
      }
    });
    mySliderTransparency.addChangeListener( new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        if( !mySliderTransparency.getValueIsAdjusting() ) {
          apply();
        }
      }
    });
    myTextHdrX.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    myTextHdrY.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    myTextMaxDistance.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
    myTextLineSize.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        apply();
      }
    });
  }
  public void addListener( csIWellPathDialogListener listener ) {
    myListeners.add( listener );
  }
  public void removeListener( csIWellPathDialogListener listener ) {
    myListeners.remove( listener );
  }
  private void addBundle() {
    WellPathBundle bundle = new WellPathBundle();
    myWellPathBundleList.add(bundle);
    updatePanel();
  }
  private void removeBundle( WellPathBundle bundle ) {
    myWellPathBundleList.remove( bundle );
    updatePanel();
  }
  public int getNumWellPaths() {
    return myWellPathBundleList.size();
  }

  private void updatePanel() {
    myMainPanel.removeAll();
    for( int ibundle = 0; ibundle < myWellPathBundleList.size(); ibundle++ ) {
      WellPathBundle bundle = myWellPathBundleList.get(ibundle);
      int xp = 0;
      myMainPanel.add(bundle.boxShow, new GridBagConstraints(
          xp++,ibundle,1,1,0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,5,0,0),0,0));
      myMainPanel.add(bundle.colorButton, new GridBagConstraints(
          xp++,ibundle,1,1,0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,5,0,0),0,0));
      myMainPanel.add(bundle.textFilename, new GridBagConstraints(
          xp++,ibundle,1,1,1.0,0,GridBagConstraints.CENTER,GridBagConstraints.HORIZONTAL,new Insets(5,10,0,0),0,0));
      myMainPanel.add(bundle.buttonOpen, new GridBagConstraints(
          xp++,ibundle,1,1,0.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,10,0,0),0,0));
      myMainPanel.add(bundle.buttonRemove, new GridBagConstraints(
          xp++,ibundle,1,1,0.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,10,0,10),0,0));
    }
    int yp = myWellPathBundleList.size();
    myMainPanel.add( myButtonAdd, new GridBagConstraints(
        0,yp++,3,1,1.0,0,GridBagConstraints.WEST,GridBagConstraints.NONE,new Insets(5,10,5,10),0,0));
    myMainPanel.add(Box.createVerticalGlue(), new GridBagConstraints(
        0,yp++,1,1,0,1.0,GridBagConstraints.WEST,GridBagConstraints.BOTH,new Insets(0,0,0,0),0,0));
    myMainPanel.revalidate();
    pack();
    repaint();
  }
  //------------------------------------------------------------------------------------------------
  //
  protected void apply() {
    apply( false );
  }
  private void apply( boolean isNewWell ) {
    csWellPathOverlayAttr overlayAttr = new csWellPathOverlayAttr();
    overlayAttr.transparency = Math.min( Math.max( mySliderTransparency.getValue(), csWellPathOverlayAttr.MIN_TRANSPARENCY ), csWellPathOverlayAttr.MAX_TRANSPARENCY );
    overlayAttr.showShadow   = myBoxShowShadow.isSelected();
    overlayAttr.shadowDepth  = Math.min( Math.max( mySliderShadowDepth.getValue(), csWellPathOverlayAttr.MIN_SHADOW_DEPTH), csWellPathOverlayAttr.MAX_SHADOW_DEPTH );
    try {
      overlayAttr.maxDisplayDistance = Integer.parseInt( myTextMaxDistance.getText() );
    }
    catch( NumberFormatException e ) {
      myTextMaxDistance.setText( overlayAttr.maxDisplayDistance + "" );
    }
    try {
      overlayAttr.lineSize    = Integer.parseInt( myTextLineSize.getText() );
    }
    catch( NumberFormatException e ) {
      myTextLineSize.setText( overlayAttr.lineSize + "" );
    }
    overlayAttr.hdrName_binx = myTextHdrX.getText();
    overlayAttr.hdrName_biny = myTextHdrY.getText();
    ArrayList<csWellPathAttr> wellPathAttrList = new ArrayList<csWellPathAttr>(myWellPathBundleList.size());
    for( WellPathBundle bundle : myWellPathBundleList ) {
      if( !bundle.attr.posXYZList.isEmpty() ) wellPathAttrList.add( bundle.attr );
    }
    fireWellPathEvent( wellPathAttrList, overlayAttr );
  }
  
  public static boolean readFile( JFrame window, csWellPathAttr attr ) {
    attr.posXYZList.clear();
    BufferedReader reader = null;
    try {
      File file = new File(attr.path);
      reader = new BufferedReader( new FileReader(file) );
    }
    catch( IOException e ) {
      if( window != null ) {
        JOptionPane.showMessageDialog( window,
            "Error occurred when trying to read in well path file\n" + attr.filename + "\n System message:\n" + e.getMessage(),
            "Error", JOptionPane.ERROR_MESSAGE );
      }
      return false;
    }
    String line;
    try {
      while( (line = reader.readLine()) != null ) {
        String[] tokens = line.replaceAll("\t"," ").trim().split(" ");
        if( tokens.length < 3 ) {
          continue; // Empty line or not enough columns
        }
        if( tokens[0].charAt(0) == '#' ) continue; // Comment line 
        try {
          csPoint3D pos = new csPoint3D();
          pos.x = Double.parseDouble(tokens[0]);
          pos.y = Double.parseDouble(tokens[1]);
          pos.z = Double.parseDouble(tokens[2]);
          attr.posXYZList.add(pos);
        } catch (NumberFormatException e) {
          e.printStackTrace();
        }
      }
      reader.close();
    }
    catch( IOException e ) {
      return false;
    }
    if( attr.posXYZList.size() == 0 ) {
      JOptionPane.showMessageDialog( window,
          "Well path file\n" + attr.filename + "\ndoes not contain any valid points\n",
          "Error", JOptionPane.ERROR_MESSAGE );
      return false;
    }
    return true;
  }
  private boolean openBundle( WellPathBundle bundle ) {
    JFileChooser fc = new JFileChooser(myCurrentDirectory);
    int option = fc.showOpenDialog( this );
    myCurrentDirectory = fc.getCurrentDirectory().getAbsolutePath();
    if( option == JFileChooser.APPROVE_OPTION ) {
      File file = fc.getSelectedFile();
      bundle.attr.path     = file.getAbsolutePath();
      bundle.attr.filename = file.getName();
      bundle.textFilename.setText( bundle.attr.filename );
      bundle.textFilename.setToolTipText( bundle.attr.path );
      bundle.attr.posXYZList.clear();
      boolean success = csWellPathDialog.readFile( myParentFrame, bundle.attr );
      if( !success ) return false;
      bundle.boxShow.setSelected(true);
      bundle.attr.showWell = true;
      apply( true );
      fireWellPathDirectoryEvent( myCurrentDirectory);
      return true;
    }
    return false;
  }
  protected void fireWellPathEvent( ArrayList<csWellPathAttr> wellPathAttrList, csWellPathOverlayAttr wellPathOverlayAttr ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).updateWellPaths( wellPathAttrList, wellPathOverlayAttr );
    }
  }
  protected void fireWellPathDirectoryEvent( String directory ) {
    for( int i = 0; i < myListeners.size(); i++ ) {
      myListeners.get(i).updateWellPathDirectory( directory );
    }
  }

  class WellPathBundle {
    JCheckBox boxShow;
    JTextField textFilename;
    JButton buttonOpen;
    JButton buttonRemove;
    csWellPathAttr attr;
    csColorButton colorButton;
    final int bundleID;
    
    WellPathBundle() {
      this( null );
    }
    WellPathBundle( String filename ) {
      bundleID = csWellPathDialog.BUNDLE_ID_COUNTER++;
      attr     = new csWellPathAttr();

      boxShow       = new JCheckBox( "", false );
      textFilename  = new JTextField( filename != null ? filename : " " );
      buttonOpen    = new JButton("Open");
      buttonRemove  = new JButton("Remove");

      textFilename.setToolTipText( filename != null ? filename : " " );
      boxShow.setToolTipText( "Show well path" );
      buttonRemove.setToolTipText( "Remove well path" );

      boxShow.setEnabled( false );
      
      Dimension dim = buttonOpen.getPreferredSize();
      textFilename.setPreferredSize( new Dimension(3*dim.width, textFilename.getPreferredSize().height) );
      textFilename.setEditable(false);

      Color color = COLOR_LIST[ bundleID % COLOR_LIST.length ];
      
      colorButton   = new csColorButton( csWellPathDialog.this, color );
      colorButton.addColorChangeListener( new csColorChangeListener() {
        public void colorChanged( Object obj, Color color ) {
//          int item = ((csColorButton)obj).getID();
          attr.color = colorButton.getColor();
          apply();
        }
      });
      
      buttonRemove.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          removeBundle( WellPathBundle.this );
          apply();
        }
      });
      boxShow.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          attr.showWell = boxShow.isSelected();
          apply();
        }
      });
      buttonOpen.addActionListener( new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
          if( openBundle( WellPathBundle.this ) ) {
            buttonOpen.setEnabled( false );
            boxShow.setEnabled( true );
          }
        }
      });
    }
  }
  
  public static void main( String[] args ) {
    JFrame frame = new JFrame();
//    JPanel panel = new JPanel();
    csWellPathDialog d = new csWellPathDialog(frame,"/users/bjolofs/tmp");
    d.setVisible(true);
  }
}
