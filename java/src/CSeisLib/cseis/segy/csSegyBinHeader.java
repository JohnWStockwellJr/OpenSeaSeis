/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.segy;

/**
 * SEGY binary header
 * @author Bjorn Olofsson
 */
public class csSegyBinHeader {
  public static final int HDR_SIZE = 400;

  public int jobID;
  public int lineNum;
  public int reelNum;
  public short numTraces;
  public short numAuxTraces;
  public short sampleIntUS;
  public short sampleIntOrigUS;
  public short numSamples;
  public short numSamplesOrig;

  /// SEGY data format
  /// 1 = 32-bit IBM floating point
  /// 2 = 32-bit fixed-point (integer)
  /// 3 = 16-bit fixed-point (integer)
  /// 4 = 32-bit fixed-point with gain code (obsolete)
  /// 5 = 32-bit IEEE floating point
  public short dataSampleFormat;
  public short fold;
  /**
  * Sorting code
  * -1 = Other
  * 0 = Unknown
  * 1 = As recorded (no sorting)
  * 2 = CDP ensemble
  * 3 = Single fold continuous profile
  * 4 = horizontally stacked
  * 5 = Common source point
  * 6 = Common receiver point
  * 7 = Common offset point
  * 8 = Common mid-point
  * 9 = Common conversion point
  */
  public short sortCode;
  /// Vertical sum code. 1 = no sum, 2 = two sum ... N=M-1 (M=2-32767)
  public short vertSumCode;
  public short sweepFreqStart; // Hz
  public short sweepFreqEnd;   // Hz
  public short sweepCode;
  public short taperType;
  public short correlatedTraces;
  public short gainRecovered;
  /// Amplitude recovery method. 1 = none, 2 = spherical divergence, 3 = AGC, 4 = other
  public short ampRecoveryMethod;
  /// Measurement system. 1 = Meters, 2 = Feet
  public short unitSystem;
  /// Impule signal polarity. 1 = Increase in pressure is negative number, 2 = otherwise
  public short polarity;
  public short vibPolarityCode;
  /// SEGY revision number. 0100 means version 1.00.
  public short revisionNum;
  public short fixedTraceLengthFlag;
  public short numExtendedBlocks;

  public csSegyBinHeader() {
  }
  public String headerDump() {
    String out = "";
    out += "Job identification number     : " + jobID + "\n";
    out += "Line number                   : " + lineNum + "\n";
    out += "Reel number                   : " + reelNum  + "\n";
    out += "Number of traces/ensemble     : " + numTraces  + "\n";
    out += "Number of aux traces/ensemble : " + numAuxTraces  + "\n";
    out += "Sample interval [us]          : " + sampleIntUS  + "\n";
    out += "Sample interval(field tape)[us]:" + sampleIntOrigUS  + "\n";
    out += "Number of samples             : " + numSamples  + "\n";
    out += "Number of samples (field tape): " + numSamplesOrig  + "\n";
    out += "Data sample format            : " + dataSampleFormat  + "\n";
    out += "Ensemble fold                 : " + fold  + "\n";
    out += "Trace sort code               : " + sortCode  + "\n";
    out += "Vertical sum code             : " + vertSumCode  + "\n";
    out += "Sweep freq at start           : " + sweepFreqStart  + "\n";
    out += "Sweep freq at end             : " + sweepFreqEnd  + "\n";
    out += "Sweep type code               : " + sweepCode  + "\n";
    out += "Taper type                    : " + taperType  + "\n";
    out += "Correlated data traces        : " + correlatedTraces  + "\n";
    out += "Binary gain recovered         : " + gainRecovered  + "\n";
    out += "Amplitude recovery method     : " + ampRecoveryMethod  + "\n";
    out += "Measurement system            : " + unitSystem  + "\n";
    out += "Impulse signal polarity       : " + polarity  + "\n";
    out += "Vibratory polarity code       : " + vibPolarityCode  + "\n";
    out += "SEG Y Format Revision number  : " + revisionNum  + "\n";
    out += "Fixed length trace flag       : " + fixedTraceLengthFlag  + "\n";
    out += "Extended header blocks        : " + numExtendedBlocks  + "\n";
    return out;
  }
}


