/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.cmapgen;

import cseis.general.csCustomColorMap;
import cseis.general.csStandard;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import javax.swing.*;

/**
 * Dialog window providing functionality to save all custom color maps as ASCII files.
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csColorMapSaveDialog extends JDialog {
  private JTextField myTextDirectory;
  private JButton myButtonSelect;
  private JCheckBox myBoxDeleteOtherFiles;
  private JButton myButtonShowOtherFiles;

  private JButton myButtonApply;
  private JButton myButtonClose;
  private String myOutputDirectory;
  private final java.util.List<csCustomColorMap> myColorMapList;

  public csColorMapSaveDialog( JFrame parentFrame, String directory, java.util.List<csCustomColorMap> cmapList ) {
    super( parentFrame, "Save all custom color maps", true );
    myOutputDirectory = directory;
    myColorMapList = cmapList;

    myBoxDeleteOtherFiles = new JCheckBox("Remove obsolete color map files");
    myBoxDeleteOtherFiles.setSelected(false);
    myTextDirectory = new JTextField( myOutputDirectory );
    myButtonSelect = new JButton("Select");
    myButtonShowOtherFiles = new JButton("Show");
    myButtonSelect.setToolTipText("Select directory");
    myTextDirectory.setToolTipText("Directory name where color maps shall be saved");
    myBoxDeleteOtherFiles.setToolTipText("Select to automatically remove all obsolete color map files from given directory");
    myButtonShowOtherFiles.setToolTipText("Show listing of obsolete color map files in given directory");
    
    JPanel panelSelect = new JPanel(new GridBagLayout());
    panelSelect.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Select output directory"),
        csStandard.INNER_EMPTY_BORDER ) );
    JPanel panelRemove = new JPanel(new GridBagLayout());
    panelRemove.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Obsolete color map files"),
        csStandard.INNER_EMPTY_BORDER ) );
    JPanel panelHelp = new JPanel(new BorderLayout());
    panelHelp.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("Info"),
        csStandard.INNER_EMPTY_BORDER ) );
    
    JTextArea area = new JTextArea();
    area.setEditable(false);
    area.setWrapStyleWord( true );
    area.setLineWrap( true );
    area.append("All custom color maps will be saved as separate ASCII files to the specified directory.\n");
    area.append("WARNING:\nAny existing color map files with the same name will be overwritten.");
    area.append("The file name of each color map is auto-generated as follows: <color map name>.cmap");
    panelHelp.add(area);
    
    int xp = 0;
    int yp = 0;
    panelSelect.add( myTextDirectory, new GridBagConstraints(
        xp++, yp, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 10, 0, 0 ), 0, 0 ) );
    panelSelect.add( myButtonSelect, new GridBagConstraints(
        xp++, yp, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 10 ), 0, 0 ) );

    panelRemove.add( myBoxDeleteOtherFiles, new GridBagConstraints(
        0, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 10, 0, 10 ), 0, 0 ) );
    panelRemove.add( myButtonShowOtherFiles, new GridBagConstraints(
        1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER,
        GridBagConstraints.NONE, new Insets( 0, 10, 0, 10 ), 0, 0 ) );
    
    myButtonApply = new JButton("Save");
    myButtonClose = new JButton("Cancel");
    this.getRootPane().setDefaultButton(myButtonApply);
    
    myButtonApply.setToolTipText("Save color maps to specified directory");
    myButtonClose.setToolTipText("Cancel operation, do not save color maps");
    
    JPanel panelButtons = new JPanel( new GridBagLayout() );
    xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonApply, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonClose, new GridBagConstraints(
        xp++, 0, 1, 1, 0.1, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.BOTH, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.35, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );

    JPanel panelTop = new JPanel(new BorderLayout());
    panelTop.add(panelSelect,BorderLayout.NORTH);
    panelTop.add(panelRemove,BorderLayout.SOUTH);
    
    yp = 0;

    panelAll.add(panelTop,BorderLayout.NORTH);
    panelAll.add(panelHelp,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);
    getContentPane().add(panelAll);

//     updateSettings( myParams, traceHeaders );
    pack();
    setSize( new Dimension(400,300) );
    setLocationRelativeTo( parentFrame );

    myButtonSelect.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        selectDirectory();
      }
    });
    myButtonApply.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        if( apply() ) {
          cancel();
        }
      }
    });
    myButtonClose.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        cancel();
      }
    });
    myButtonShowOtherFiles.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        showObsoleteColorMapFiles();
      }
    });
    myTextDirectory.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        myOutputDirectory = myTextDirectory.getText();
      }
    });

    this.addWindowListener(new WindowAdapter() {
      @Override
      public void windowClosing(WindowEvent e) {
        cancel();
      }
    });
  }
  private void selectDirectory() {
    JFileChooser fc = new JFileChooser( myOutputDirectory );
    fc.setFileSelectionMode( JFileChooser.DIRECTORIES_ONLY );
    int option = fc.showOpenDialog( this );
    if( option == JFileChooser.APPROVE_OPTION ) {
      myOutputDirectory = fc.getSelectedFile().getAbsolutePath();
      myTextDirectory.setText(myOutputDirectory);
    }
  }
  private void cancel() {
    dispose();
  }
  private boolean apply() {
    String directoryName = myTextDirectory.getText();
    boolean deleteOtherFiles = myBoxDeleteOtherFiles.isSelected();
    String message = "";
    if( deleteOtherFiles ) {
      java.util.List<String> filenameList = retrieveObsoleteColorMapFiles(directoryName);
      int counter = 0;
      for( int i = 0; i < filenameList.size(); i++ ) {
        File file = new File( directoryName + java.io.File.separatorChar + filenameList.get(i) );
        boolean success = file.delete();
        if( success ) counter += 1;
      }
      if( counter > 0 ) message += "Deleted " + counter + " obsolete color map files.\n";
    }
    String filename = "";
    try {
      for( int imap = 0; imap < myColorMapList.size(); imap++ ) {
        csCustomColorMap cmap = myColorMapList.get(imap);
        filename = directoryName + java.io.File.separatorChar + cmap.toString() + ".cmap";
        csColorMapGenerator.saveColorMap( cmap, filename );
      }
    }
    catch( FileNotFoundException ex ) {
      JOptionPane.showMessageDialog(
              this,
              "File " + filename + ": \n" + ex.getMessage(),
              "Error", JOptionPane.ERROR_MESSAGE);
      return false;
    }
    catch( IOException ex ) {
      JOptionPane.showMessageDialog(
              this,
              "File " + filename + ": \n" + ex.getMessage(),
              "Error", JOptionPane.ERROR_MESSAGE);
      return false;
    }
    JOptionPane.showMessageDialog(
            this,
            message +
            "Saved " + myColorMapList.size() + " color maps to directory...\n  " + directoryName,
            "Info", JOptionPane.INFORMATION_MESSAGE);
    return true;
  }
  public java.util.List<String> retrieveObsoleteColorMapFiles( String directoryName ) {
    ArrayList<String> newFilenameList = new ArrayList<String>();
    for( int imap = 0; imap < myColorMapList.size(); imap++ ) {
      csCustomColorMap cmap = myColorMapList.get(imap);
      newFilenameList.add( cmap.toString() + ".cmap" );
    }
    File dirFile = new File( directoryName );
    File[] fileList = dirFile.listFiles( new FilenameFilter() {
      @Override
      public boolean accept( File dir, String name ) {
        return name.endsWith(".cmap");
      }
    });
    ArrayList<String> obsoleteFilenameList = new ArrayList<String>();
    for( int i = 0; i < fileList.length; i++ ) {
      String filename = fileList[i].getName();
      if( !newFilenameList.contains(filename) ) {
        obsoleteFilenameList.add(filename);
      }
    }
    return obsoleteFilenameList;
  }
  private void showObsoleteColorMapFiles() {
    String directoryName  = myTextDirectory.getText();
    java.util.List<String> filenameList = retrieveObsoleteColorMapFiles(directoryName);
    String message;
    if( filenameList.isEmpty() ) {
      message = "No obsolete color map files found in directory'" + directoryName + "'";
    }
    else {
      message = "List of obsolete color map files in directory '" + directoryName + "':\n";
      for( int i = 0; i < filenameList.size(); i++ ) {
        String filename = filenameList.get(i);
        message += "   " + filename + "\n";
      }
    }
    JOptionPane.showMessageDialog( this, message, "Onsolete color map files", JOptionPane.INFORMATION_MESSAGE);
  }
}


