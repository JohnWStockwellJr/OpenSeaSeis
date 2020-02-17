/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.velanal;

import cseis.general.csColorMap;
import cseis.seaview.SeaView;
import cseis.seaview.csISeisPaneBundleMouseModeListener;
import cseis.seaview.csSeisPaneBundle;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seisdisp.csAnnotationAttributes;
import cseis.seisdisp.csISeisViewListener;
import cseis.seisdisp.csMouseModes;
import cseis.seisdisp.csSeisDispSettings;
import cseis.velocity.csEnsemble;
import cseis.velocity.csVelEnsembleInfo;
import cseis.velocity.csVelFunction;
import cseis.velocity.csVelPickOverlay;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

public class csVelocityAnalysisManager implements csIVelSetupDialogListener, csIVelPickDialogListener {
  private final SeaView myParentFrame;
  private csSeisPaneBundle mySemblanceBundle;
  private csVelocityAnalysisSetupDialog mySetupDialog;
  private csVelocityPickingDialog myVelPickDialog;
  private csVelPickOverlay myVelPickOverlay;
  private boolean myHasBeenSetup = false;
  private csVelEnsembleInfo myEnsInfo;
  private String myVelFieldDirPath;
  
  private JFileChooser myFileChooser;
  
  public csVelocityAnalysisManager( SeaView parentFrame, csSeisPaneBundle[] bundleList, csSeisPaneBundle semblanceBundleProposed, String velFieldDirPath ) {
    myParentFrame = parentFrame;
    myHasBeenSetup = false;
    mySetupDialog = new csVelocityAnalysisSetupDialog( parentFrame, bundleList, semblanceBundleProposed, this );
    myVelFieldDirPath = velFieldDirPath;
  }
  public File getVelFieldDir() {
    if( myFileChooser != null ) {
      return myFileChooser.getCurrentDirectory();
    }
    else if( myVelFieldDirPath == null ) {
      return null;
    }
    else {
      return( new File(myVelFieldDirPath) );
    }
  }
  public void open() {
    if( !myHasBeenSetup ) {
      mySetupDialog.setVisible(true);
    }
    if( myHasBeenSetup && mySemblanceBundle != null ) {
      if( myVelPickDialog == null ) {
        myVelPickDialog = new csVelocityPickingDialog( myParentFrame, myVelPickOverlay, mySemblanceBundle.title, this );
        myVelPickDialog.refresh();
      }
      myVelPickDialog.setVisible(true);
      myVelPickOverlay.activate();
      myVelPickOverlay.refresh();
      mySemblanceBundle.setLocalMouseMode( csMouseModes.VELPICK_MODE );
    }
  }
  public void stop() {
    if( myVelPickDialog != null ) {
      myVelPickDialog.dispose();
    }
    velPickDialogClosing();
  }
  @Override
  public void setupVelocityAnalysis( csSeisPaneBundle bundle, csVelEnsembleInfo ensInfo ) {
    myEnsInfo = ensInfo;
    try {
      if( bundle == null || ensInfo == null ) throw( new Exception("Unknown error occurred (-->seismic bundle and/or ensemble info are null)") );
      mySemblanceBundle = bundle;
      csHeaderDef[] hdef = mySemblanceBundle.getHeaderDef();
      for( int ihdr = 0; ihdr < hdef.length; ihdr++ ) {
        String name = hdef[ihdr].name;
        if( name.compareTo( ensInfo.hdrNameVel ) == 0 ) {
          ensInfo.hdrIndexVel = ihdr;
          if( ensInfo.getDim() == csEnsemble.DIM_1D ) break;
        }
        else if( ensInfo.getDim() != csEnsemble.DIM_1D && name.compareTo( ensInfo.hdrNameEns1 ) == 0 ) {
          ensInfo.hdrIndexEns1 = ihdr;
        }
        else if( ensInfo.getDim() == csEnsemble.DIM_3D && name.compareTo( ensInfo.hdrNameEns2 ) == 0 ) {
          ensInfo.hdrIndexEns2 = ihdr;
        }
      }
      myVelPickOverlay = new csVelPickOverlay( mySemblanceBundle.seisView, ensInfo );
      myVelPickOverlay.traceBufferChanged( mySemblanceBundle.seisView.getTraceBuffer() );
      mySemblanceBundle.setLocalMouseMode( csMouseModes.VELPICK_MODE );
      mySemblanceBundle.seisView.addSeisViewListener( new csISeisViewListener() {
        @Override
        public void traceBufferChanged( csISeismicTraceBuffer traceBuffer ) {
          if( myVelPickOverlay != null ) {
            try {
              myVelPickOverlay.traceBufferChanged( traceBuffer );
            }
            catch( Exception exc ) {
              JOptionPane.showMessageDialog( myParentFrame, "Error occurred when refreshing semblance for velocity picking:\n" +
              exc.getMessage() , "Error", JOptionPane.ERROR_MESSAGE );
            }
          }
        }
        @Override
        public void changedSettings( csSeisDispSettings settings ) {
        }
        @Override
        public void vertScrollChanged( int scrollValue ) {
        }
        @Override
        public void horzScrollChanged( int scrollValue ) {
        }
        @Override
        public void sizeChanged( Dimension size ) {
        }
        @Override
        public void mouseExited() {
        }
      });
      mySemblanceBundle.addSeisPaneBundleMouseModeListener( new csISeisPaneBundleMouseModeListener() {
        @Override
        public void updateBundleMouseMode(int mouseMode) {
          if( mouseMode != csMouseModes.VELPICK_MODE ) {
            if( myVelPickDialog != null && myVelPickDialog.isVisible() ) myVelPickDialog.dispose();
            if( myVelPickOverlay != null ) {
              myVelPickOverlay.deactivate();
              myVelPickOverlay.refresh();
            }
          }
        }
      });
      //------------------------------
      // Apply default display settings to semblance display
      mySemblanceBundle.seisView.addOverlay( myVelPickOverlay );
      csSeisDispSettings settings = mySemblanceBundle.seisView.getDispSettings();
      settings.isVIDisplay = true;
      settings.showWiggle  = false;
      settings.isPosFill   = false;
      settings.isNegFill   = false;
      settings.viType      = csSeisDispSettings.VA_TYPE_2DSPLINE;
      settings.viColorMap.resetColors( new csColorMap(csColorMap.COLD_WARM,settings.viColorMap.getColorMapType()) );
      settings.scaleType   = csSeisDispSettings.SCALE_TYPE_RANGE;
      mySemblanceBundle.seisView.updateDispSettings( settings );
      
      csAnnotationAttributes ann_attr = new csAnnotationAttributes();
      ann_attr.omitRepeating = false;
      ann_attr.showSequential = false;
      ann_attr.fixedTraceLabelStep = false;
      ann_attr.traceLabelStep = 60;
      int numHeaders = 1;
      if( ensInfo.getDim() != csEnsemble.DIM_1D ) numHeaders += 1;
      if( ensInfo.getDim() == csEnsemble.DIM_3D ) numHeaders += 1;
      csHeaderDef[] headers = new csHeaderDef[ numHeaders ];
      headers[0] = new csHeaderDef( ensInfo.hdrNameVel, "Velocity", csHeaderDef.TYPE_FLOAT );
      int counter = 1;
      if( ensInfo.getDim() == csEnsemble.DIM_3D ) headers[counter++] = new csHeaderDef( ensInfo.hdrNameEns2, "Ensemble2", csHeaderDef.TYPE_INT );
      if( ensInfo.getDim() != csEnsemble.DIM_1D ) headers[counter] = new csHeaderDef( ensInfo.hdrNameEns1, "Ensemble1", csHeaderDef.TYPE_INT );
      mySemblanceBundle.seisPane.updateAnnotationTraceHeaders( headers, ann_attr ); 
      mySemblanceBundle.seisPane.repaint();
      
      myHasBeenSetup = true;
    }
    catch( Exception exc ) {
      JOptionPane.showMessageDialog( mySetupDialog, "Error occurred when trying to set up semblance velocity picking:\n" +
        exc.getMessage() , "Error", JOptionPane.ERROR_MESSAGE );
      myVelPickOverlay  = null;
      mySemblanceBundle = null;
    }
  }
  private void createFileChooser() {
    myFileChooser = new JFileChooser(  );
    if( myVelFieldDirPath != null ) myFileChooser.setCurrentDirectory( new File(myVelFieldDirPath) );
  }
  @Override
  public void loadPicks() {
    csVelEnsembleInfo ensInfo = myVelPickOverlay.getEnsembleInfo();
    if( myFileChooser == null ) createFileChooser();
    int option = myFileChooser.showOpenDialog( myVelPickDialog );
    if( option == JFileChooser.APPROVE_OPTION ) {
      File file = myFileChooser.getSelectedFile();
      try {
        BufferedReader reader = new BufferedReader( new FileReader(file) );
        String line;
        ArrayList<Float> timeList = new ArrayList<Float>();
        ArrayList<Float> velList  = new ArrayList<Float>();
        ArrayList<Integer> ens1List = new ArrayList<Integer>();
        ArrayList<Integer> ens2List = new ArrayList<Integer>();
        while( (line = reader.readLine()) != null ) {
          String[] tokens = line.replaceAll("\t"," ").trim().split(" ");
          if( tokens.length > 0 && tokens[0].charAt(0) == '#' ) continue; // Comment line 
          if( tokens.length >= 2 ) {
            int index = 0;
            if( ensInfo.getDim() == csEnsemble.DIM_2D ) {
              if( tokens.length < 3 ) throw( new IOException("Too few values/columns in input file: Need 'ensemble', 'time' and 'vel'") );
              ens1List.add( (int)Double.parseDouble( tokens[index++] ) );
            }
            else if( ensInfo.getDim() == csEnsemble.DIM_3D ) {
              if( tokens.length < 4 ) throw( new IOException("Too few values/columns in input file: Need 'ensemble1', 'ensemble2', 'time' and 'vel'") );
              ens1List.add( (int)Double.parseDouble( tokens[index++] ) );
              ens2List.add( (int)Double.parseDouble( tokens[index++] ) );
            }
            timeList.add( Float.parseFloat( tokens[index++] ) );
            velList.add( Float.parseFloat( tokens[index++] ) );
          }
        }
        reader.close();
        if( timeList.isEmpty() ) return;
//-------------------------------------------------------------------
// START create vel function list
        ArrayList<csVelFunction> velFuncList = new ArrayList<csVelFunction>();
        int counter = 0;
        csEnsemble ensCurrent;
        if( ensInfo.getDim() == csEnsemble.DIM_1D ) {
          ensCurrent = new csEnsemble();
        }
        else if( ensInfo.getDim() == csEnsemble.DIM_2D ) {
          ensCurrent = new csEnsemble( ens1List.get(counter) );
        }
        else { //if( ensInfo.getDim() == csEnsemble.DIM_3D ) {
          ensCurrent = new csEnsemble( ens1List.get(counter), ens2List.get(counter) );
        }
        csVelFunction velFuncCurrent = new csVelFunction(ensCurrent);
        while( counter < timeList.size() ) {
          csEnsemble ens;
          if( ensInfo.getDim() == csEnsemble.DIM_1D ) {
            ens = new csEnsemble();
          }
          else if( ensInfo.getDim() == csEnsemble.DIM_2D ) {
            ens = new csEnsemble( ens1List.get(counter) );
          }
          else { //if( ensInfo.getDim() == csEnsemble.DIM_3D ) {
            ens = new csEnsemble( ens1List.get(counter), ens2List.get(counter) );
          }
          if( !ens.equals(ensCurrent) ) {
            velFuncList.add(velFuncCurrent);
            ensCurrent = ens;
            velFuncCurrent = new csVelFunction(ensCurrent);
          }
          velFuncCurrent.addPick(timeList.get(counter), velList.get(counter));
          counter += 1;
        }
        velFuncList.add(velFuncCurrent);
// END create vel function list
//-----------------------------------------------------------
        String name = file.getAbsolutePath();
        int separatorIndex = name.lastIndexOf( java.io.File.separatorChar );
        if( separatorIndex < 0 ) separatorIndex = 0;
        int pointIndex = name.lastIndexOf('.');
        if( pointIndex >= 0 ) name = name.substring( separatorIndex+1, pointIndex );
        myVelPickDialog.load( velFuncList, name );
        String ensText;
        if( ensInfo.getDim() == csEnsemble.DIM_1D ) {
          ensText = "N/A";
        }
        else {
          ensText = ensInfo.hdrNameEns1 + " " + ens1List.get(0) + "-" + ens1List.get(ens1List.size()-1);
          if( ensInfo.getDim() == csEnsemble.DIM_3D ) {
            ensText += " / " + ensInfo.hdrNameEns2 + " " + ens2List.get(0) + "-" + ens2List.get(ens1List.size()-1);
          }
        }
        JOptionPane.showMessageDialog( myVelPickDialog,
          "Loaded velocity field/picks from file\n'" + file.getAbsolutePath() +
            "'\nEnsemble range in file: " + ensText,
          "Info", JOptionPane.INFORMATION_MESSAGE);
      }
      catch (IOException ex) {
        JOptionPane.showMessageDialog( myVelPickDialog, "Error occurred when loading time picks:\n" + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
      }
      catch (NumberFormatException ex) {
        JOptionPane.showMessageDialog( myVelPickDialog, "Error occurred when loading time picks:\n" + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
      }
      
    }
  }

  @Override
  public void velPickDialogClosing() {
    if( mySemblanceBundle != null ) mySemblanceBundle.setLocalMouseMode( csMouseModes.NO_MODE );
  }
  @Override
  public void velPickDialogReset() {
    myHasBeenSetup = false;
    mySemblanceBundle.seisView.removeOverlay( myVelPickOverlay );
    mySemblanceBundle.setLocalMouseMode( csMouseModes.NO_MODE );
    javax.swing.SwingUtilities.invokeLater( new Runnable() {
      @Override
      public void run() {
        if( myVelPickDialog != null ) myVelPickDialog.dispose();
        myVelPickDialog = null;
        myParentFrame.openVelocityAnalysis( true );
      }
    });
  }
  @Override
  public void savePicks( java.util.Collection<csVelFunction> velFuncList, String name ) {
    if( velFuncList.isEmpty() ) {
      JOptionPane.showMessageDialog(myParentFrame, "Cannot save velocity picks.\nVelocity field " + name + " is empty.", "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }
    if( myFileChooser == null ) createFileChooser();
    myFileChooser.setSelectedFile( new File(name + ".table") );
    int option = myFileChooser.showSaveDialog( myVelPickDialog );
    if( option == JFileChooser.APPROVE_OPTION ) {
      File file = myFileChooser.getSelectedFile();
      try {
        BufferedWriter writer = new BufferedWriter( new FileWriter(file) );
        DecimalFormat floatFormat = new DecimalFormat("0.0");
        Iterator<csVelFunction> iterator = velFuncList.iterator();
        boolean isFirst = true;
        while( iterator.hasNext() ) {
          csVelFunction velFunc = iterator.next();
          if( isFirst ) {
            csEnsemble ens = velFunc.getEns();
            String text = "# SeaView " + ens.dimToString() + " velocity function \n";
            text += "# Columns: ";
            if( !velFunc.getEns().is1D() ) text += myEnsInfo.hdrNameEns1;
            if( velFunc.getEns().is3D() ) text += " " + myEnsInfo.hdrNameEns2;
            text += " time[ms] velocity[m/s]\n";
            writer.write( text );
            isFirst = false;
          }
          String ensText = velFunc.getEns().toString();
          String text = "";
          java.util.List<csVelFunction.Item> pickItemList = velFunc.getPicks();
          for( csVelFunction.Item item : pickItemList ) {
            text += ensText + " " + floatFormat.format(item.time) + " " + floatFormat.format(item.vel) + "\n";
          }
          writer.write( text );
        }
        writer.close();
        JOptionPane.showMessageDialog(myParentFrame, "Saved velocity picks to file\n" + file.getAbsolutePath(), "Info", JOptionPane.INFORMATION_MESSAGE);
      }
      catch (IOException ex) {
        JOptionPane.showMessageDialog(myParentFrame, "Error occurred when writing velocity picks:\n" + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
      }
    }
  }

}
