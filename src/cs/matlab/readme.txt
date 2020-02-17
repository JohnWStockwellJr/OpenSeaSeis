MATLAB Deployment Checklist
===========================

Use this checklist to successfully deploy your application, component, 
or library.

You can distribute a MATLAB Compiler-generated standalone application, 
component, or library to any target machine with the same operating 
system as the machine on which the application was compiled 
(the "source" development machine).

For more information, see the documentation for the MATLAB Compiler at:
http://www.mathworks.com/access/helpdesk/help/toolbox/compiler/

NOTE: Replace the items enclosed in angle brackets ( < > ) with your own 
values.

1. Ensure that the MATLAB Compiler Runtime (MCR) is installed on target 
   machines, and ensure you have installed the correct version. To verify 
   this, do the following: 

    A. Attempt to verify the MCR exists on your system.

    The MCR usually resides in these locations on these platforms:

    Operating System        Path to MCR
    ================        ===========
    Windows              C:\Program Files\MATLAB\MATLAB Compiler Runtime\v78
    Mac                  Applications/MATLAB/MATLAB_Compiler_Runtime/v78 
    Linux/Solaris        *** No predetermined installation location ***
                         Contact your system administrator.  

        NOTE: For information about the installation process and the MCR, 
              type "help MCR" or "help mcrinstaller" at the MATLAB 
              command prompt. 

    B. MATLAB 2008a uses MCR version 7.8. To verify
       the version number of the installed MCR, type the following 
       command:
       [mcrmajor,mcrminor]=mcrversion
       at the MATLAB command prompt.


2. Add the MCR directory to the path specified by the target system's 
   environment variable. 

    A. Locate the name of the environment variable to set, using the table 
       below:

    Operating System        Environment Variable
    ================        ====================
    Windows                PATH
    Linux                  LD_LIBRARY_PATH
    Solaris                LD_LIBRARY_PATH
    Mac                    DYLD_LIBRARY_PATH

    B. Set the path by doing one of the following:

        On Windows systems:

        * Add the MCR directory to the environment variable by opening a 
        command prompt and issuing the DOS command, specifying either 
        win32 or win64:

            set PATH=C:\Program Files\MATLAB\MATLAB Compiler Runtime\v78\runtime\{win32|win64};%PATH% 

        Alternately, for Windows, add the following pathname:
            C:\Program Files\MATLAB\MATLAB Compiler Runtime\v78\runtime\{win32|win64}
        to the PATH environment variable, specifying either win32 or win64, 
        by doing the following:
            1. Select the My Computer icon on your desktop.
            2. Right-click the icon and select Properties from the menu.
            3. Select the Advanced tab.
            4. Click Environment Variables.  

        On UNIX systems:

        * Add the MCR directory to the environment variable by issuing 
            the following commands:

        Linux
            setenv LD_LIBRARY_PATH
                <mcr_root>/v78/runtime/glnx86:
                <mcr_root>/v78/sys/os/glnx86:
                <mcr_root>/v78/sys/java/jre/glnx86/jre1.6.0/lib/i386/native_threads:
                <mcr_root>/v78/sys/java/jre/glnx86/jre1.6.0/lib/i386/server:
                <mcr_root>/v78/sys/java/jre/glnx86/jre1.6.0/lib/i386
            setenv XAPPLRESDIR <mcr_root>/v78/X11/app-defaults

        Solaris64
            setenv LD_LIBRARY_PATH /usr/lib/lwp:
                <mcr_root>/v78/runtime/sol64:
                <mcr_root>/v78/sys/os/sol64:
                <mcr_root>/v78/sys/java/jre/sol64/jre1.6.0/lib/sparcv9/native_threads:
                <mcr_root>/v78/sys/java/jre/sol64/jre1.6.0/lib/sparcv9/server:
                <mcr_root>/v78/sys/java/jre/sol64/jre1.6.0/lib/sparcv9
            setenv XAPPLRESDIR <mcr_root>/v78/X11/app-defaults

        Linux x86-64
            setenv LD_LIBRARY_PATH
                <mcr_root>/v78/runtime/glnxa64:
                <mcr_root>/v78/sys/os/glnxa64:
                <mcr_root>/v78/sys/java/jre/glnxa64/jre1.6.0/lib/amd64/native_threads:
                <mcr_root>/v78/sys/java/jre/glnxa64/jre1.6.0/lib/amd64/server:
                <mcr_root>/v78/sys/java/jre/glnxa64/jre1.6.0/lib/amd64 
            setenv XAPPLRESDIR <mcr_root>/v78/X11/app-defaults

        Mac
            setenv DYLD_LIBRARY_PATH
                <mcr_root>/version/runtime/maci:
                <mcr_root>/version/sys/os/maci:
                <mcr_root>/version/bin/maci:
                /System/Library/Frameworks/JavaVM.framework/JavaVM:
                /System/Library/Frameworks/JavaVM.framework/Libraries
            setenv XAPPLRESDIR <mcr_root>/version/X11/app-defaults


        NOTE: To make theses changes persistent after logout on UNIX or 
              Mac machines, modify the .cshrc file to include this setenv 
              command.
        NOTE: On Windows, the Environment Variable syntax utilizes 
              backslashes (\), delimited by semi-colons (;). On UNIX,
              the EV syntax utilizes forward slashes (/), delimited by 
              colons (:). 
        NOTE: When deploying C and C++ standalone applications, it is 
              possible to run the shell script file run_libmatlabTest.sh 
              on UNIX and Mac instead of setting environment variables. 
              See "Requirements for Standalone Applications," below.    


3. Collect the following files, based on what you are deploying: 

    Requirements for Standalone Applications:
    ============================================================

    All Platforms:
        * libmatlabTest.ctf (If present)(Component Technology File archive,
          platform-dependent on target machine)

    Windows:
        * MCRInstaller.exe (self-extracting MATLAB Compiler Runtime
          installer)
        * libmatlabTest.exe (application created by MATLAB Compiler)

    UNIX:
        * MCRInstaller.bin (MATLAB Compiler Runtime installer)
        * libmatlabTest (application created by MATLAB Compiler)
        * run_libmatlabTest.sh (shell script that can be run to 
          temporarily set environment paths and execute the application)

    Mac:
        * MCRInstaller.dmg (MATLAB Compiler Runtime installer)
        * libmatlabTest (application created by MATLAB Compiler)
        * run_libmatlabTest.sh (shell script that can be run to 
          temporarily set environment paths and execute the application)


    Requirements for C and C++ shared libraries:
    ============================================================

    All Platforms:
        * libmatlabTest.so (shared library, file extension varies by 
          platform)
        * libmatlabTest.h (library header file)
        * libmatlabTest.ctf (If present)(Component Technology File archive,
          platform-dependent on target machine, contains compiled M-code)

    Windows:
        * MCRInstaller.exe (self-extracting MATLAB Compiler Runtime
          installer)  

    UNIX:
        * MCRInstaller.bin (MATLAB Compiler Runtime installer)

    Mac:
        * MCRInstaller.dmg (MATLAB Compiler Runtime installer)


    Requirements for .NET Components:
    ============================================================

    Windows:
        * libmatlabTest.ctf (Component Technology File archive, contains 
          compiled M-code)    
        * libmatlabTest.xml (documentation files)
        * libmatlabTest.pdb (program database file - if DEBUG option 
          selected)
        * libmatlabTest.dll (component assembly file)
        * MCRInstaller.exe (self-extracting MATLAB Compiler Runtime
          archive)


    Requirements for COM Components:
    ============================================================

    Windows:
        * libmatlabTest.ctf (Component Technology File archive, contains 
          compiled M-code)
        * libmatlabTest_1_0.dll (COM component)
        * MCRInstaller.exe (self-extracting MATLAB Compiler Runtime
          archive)      


    Requirements for Java Components: 
    ============================================================

    All Platforms: 
        * libmatlabTest.jar (Java class files and compiled M-code)

    Windows:
        * MCRInstaller.exe (self-extracting MATLAB Compiler Runtime
          installer)  

    UNIX:
        * MCRInstaller.bin (MATLAB Compiler Runtime installer)

    Mac:
        * MCRInstaller.dmg (MATLAB Compiler Runtime installer)


    Requirements for COM Components for Use with Microsoft Excel:
    ============================================================

    Windows:
        * libmatlabTest_1_0.dll (COM component)
        * libmatlabTest.ctf (Component Technology File archive, contains 
          compiled M-code)
        * libmatlabTest.xla (Optionally created and placed in
           directory by the user)  
        * MCRInstaller.exe (self-extracting MATLAB Compiler Runtime
          archive)

