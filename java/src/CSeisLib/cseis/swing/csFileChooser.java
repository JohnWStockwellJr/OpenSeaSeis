/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.awt.Component;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

/**
 * Customised file chooser
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csFileChooser extends JFileChooser {
  public csFileChooser() {
  }
  public csFileChooser( String title ) {
    setDialogTitle(title);
  }
  /**
   * Select file for writing
   * @param parentComponent
   * @return true if file was successfully selected
   */
  public boolean selectFileForWriting( Component parentComponent ) {
    int option = showSaveDialog( parentComponent );
    while( option == JFileChooser.APPROVE_OPTION && getSelectedFile().exists() ) {
      String[] buttonLabels = {"Yes", "No"};
      int selectedOption = JOptionPane.showOptionDialog( parentComponent,
          "The selected file exists already.\n" +
          "'" + getSelectedFile().getName() + "'\n" +
          "OK to overwrite existing file?",
          "Confirm overwrite",
          JOptionPane.YES_NO_OPTION,
          JOptionPane.QUESTION_MESSAGE,
          null, buttonLabels, buttonLabels[1] );
      if( selectedOption == JOptionPane.YES_OPTION ) {
        break;
      }
      option = showSaveDialog( parentComponent );
    }
    return( option == JFileChooser.APPROVE_OPTION );
  }
}


