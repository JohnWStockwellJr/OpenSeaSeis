/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Class to test csDockPaneManager.
 */
@SuppressWarnings("serial")
public class TestDockManager extends JFrame {
  public TestDockManager() {
    super("Test dockable manager");
    final csDockPaneManager manager = new csDockPaneManager();
    final JComboBox<String> combo = new JComboBox<String>();
    combo.addItem( "One row" );
    combo.addItem( "One col" );
    combo.addItem( "Two rows" );
    combo.addItem( "Two cols" );
    combo.addItem( "Three rows" );
    combo.addItem( "Three cols" );
    combo.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        int index = combo.getSelectedIndex();
        System.out.println("Selected index: " + index);
        if( index >= 0 ) {
          int layout = manager.resetWindowLayout( index+1, true );
          System.out.println("Layout = " + layout);
        }
      }
    });
    int numPanes = 8;
    final JTextField text = new JTextField(""+numPanes);
    text.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        int numPanes = Integer.parseInt(text.getText());
        System.out.println("Num panes: " + numPanes);
        if( numPanes > 0 ) {
          manager.setNumPanes(numPanes);
        }
      }
    });
    JMenuItem menuView = new JMenuItem("Select panels");
    JMenu menu = new JMenu("View");
    menu.add(menuView);
    JMenuBar menuBar = new JMenuBar();
    this.setJMenuBar(menuBar);
    menuBar.add(menu);
    menuView.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        final JDialog dialog = new JDialog(TestDockManager.this,"Panel selection",true);
        JPanel panelButton = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        JButton buttonClose = new JButton("Close");
        panelButton.add( buttonClose );
        dialog.getContentPane().add( panelButton, BorderLayout.SOUTH );
        dialog.getContentPane().add( manager.createSelectionPanel(), BorderLayout.CENTER );
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.pack();
        dialog.setLocationRelativeTo(TestDockManager.this);
        buttonClose.addActionListener(new ActionListener() {
          @Override
          public void actionPerformed(ActionEvent e) {
            dialog.dispose();
          }
        });
        dialog.setVisible(true);
      }
    });
    
    Color[] colors = new Color[numPanes];
    colors[0] = Color.red;
    colors[1] = Color.green;
    colors[2] = Color.white;
    colors[3] = Color.yellow;
    colors[4] = Color.blue;
    colors[5] = Color.pink;
    colors[6] = Color.gray;
    colors[7] = Color.orange;
//    colors[3] = Color.red;numPanes
    JPanel panelTop = new JPanel(new FlowLayout());
    panelTop.add(combo);
    panelTop.add(text);
    csDockPaneButtonSelection bselect = new csDockPaneButtonSelection();
    bselect.scrollbars = false;
    for( int i = 0; i < numPanes; i++ ) {
      JPanel panel = new JPanel(new GridBagLayout());
      panel.setBorder( BorderFactory.createCompoundBorder(
          BorderFactory.createTitledBorder("Panel #" + (i+1)),
          BorderFactory.createEmptyBorder(2, 2, 2, 2) ) );
      panel.add( Box.createVerticalGlue(), new GridBagConstraints(
          0, 0, 1, 1, 1.0, 0.5, GridBagConstraints.WEST,
          GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      panel.add( new JLabel("LABEL #" + (i+1)), new GridBagConstraints(
          0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER,
          GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      panel.add( Box.createVerticalGlue(), new GridBagConstraints(
          0, 2, 1, 1, 1.0, 0.5, GridBagConstraints.WEST,
          GridBagConstraints.BOTH, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
      panel.setBackground(colors[i]);
//      manager.addPane(panel,false);
      if( i == 3 ) {
        csDockPane dockPane = new csDockPane( manager, panel, "Panel #" + (i+1) +" this_is_a_long_name for the_panel", bselect );
        manager.addDockPane(dockPane,false);
      }
      else {
        csDockPane dockPane = new csDockPane( manager, panel, "Panel #" + (i+1), bselect );
        manager.addDockPane(dockPane,false);
      }
    }
    getContentPane().add(panelTop,BorderLayout.NORTH);
    getContentPane().add(manager);
    pack();
    setSize(1000,600);
    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  }
  public static void main(String[] args) {
    TestDockManager test = new TestDockManager();
    test.setVisible(true);
  }
  
}

