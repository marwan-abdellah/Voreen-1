####################################################################################################
# CMake configuration file for Voreen.                                                             #
#                                                                                                  #
# It is usually not necessary to edit this file, since build options are cached by CMake           #
# and can be conveniently specified via the CMake GUI or ccmake. However, if you want to           #
# keep settings permanently, even across deletions of the CMake cache, you may copy                #
# this file to 'config.cmake' and edit it.                                                         #
####################################################################################################

####################################################################################################
#  Projects                                                                                        
####################################################################################################

# Main Projects (enabled by default)
OPTION(VRN_BUILD_LIB_TGT        "Build tgt library?"                                            ON )
OPTION(VRN_BUILD_LIB_VOREENCORE "Build voreen_core library?"                                    ON )
OPTION(VRN_BUILD_LIB_VOREENQT   "Build voreen_qt library?"                                      ON )
OPTION(VRN_BUILD_VOREENVE       "Build VoreenVE application?"                                   ON )
OPTION(VRN_BUILD_VOREENTOOL     "Build VoreenTool console application?"                         ON )

# Additional Projects (disabled by default)
OPTION(VRN_BUILD_SIMPLEQT       "Build Qt-based sample application?"                            OFF)
OPTION(VRN_BUILD_SIMPLEGLUT     "Build GLUT-based sample application?"                          OFF)

# Test apps
IF(EXISTS ${VRN_HOME}/apps/tests)
    OPTION(VRN_BUILD_TESTAPPS   "Build Voreen testing applications?"                            OFF)
ENDIF()
    
# Advanced Projects (mainly for internal use)
IF(EXISTS ${VRN_HOME}/apps/itk_wrapper)
    OPTION(VRN_BUILD_ITKWRAPPER     "Build ITK-wrapper application?"                            OFF)
    MARK_AS_ADVANCED(VRN_BUILD_ITKWRAPPER)
ENDIF()

####################################################################################################
#  Modules
####################################################################################################
SET(VRN_CUSTOM_MODULEDIR        ${VRN_HOME}/custommodules     
    CACHE PATH "Directory where custom modules are stored")

####################################################################################################
#  Options
####################################################################################################
SET(VRN_BINARY_OUTPUT_DIR       ${VRN_HOME}/bin
    CACHE PATH "Directory where binary files (apps and libraries) are put in" )
OPTION(VRN_SHARED_LIBS          "Build shared libraries?"                                       ON )
OPTION(VRN_PRECOMPILED_HEADER   "Use pre-compiled headers?"                                     ON )
OPTION(VRN_ADD_INSTALL_TARGET   "Add install target for creating a binary distribution?"        OFF)

IF(WIN32)
    OPTION(VRN_GROUP_SOURCE_FILES   "Group source files by sub directory?"                      ON )
    OPTION(VRN_COPY_EXTERNAL_DLLS   "Copy external DLLs to bin directory?"                      ON )
    OPTION(VRN_INCREMENTAL_LINKING  "Enable incremental linking in Visual Studio debug builds?" ON )
    OPTION(VRN_GENERATE_MANIFEST    "Generate manifest in Visual Studio debug builds?"          OFF)
    OPTION(VRN_DISABLE_ZERO_CHECK   
        "Do not generate ZERO_CHECK target (might cause unnecessary PCH re-compilations)"       ON )

    OPTION(VRN_DEPLOYMENT   "Settings/install target for deployment"                           OFF)
    MARK_AS_ADVANCED(VRN_DEPLOYMENT)
ENDIF()
