include_directories( $ENV{LARCV_INCDIR} )
set( BOOST_LIB -L$ENV{BOOST_LIB} -lboost_system )
set( LARCV_LIB -L$ENV{LARCV_LIBDIR} -lLArCVCoreDataFormat )

art_make( MODULE_LIBRARIES larcorealg_Geometry
	  larcore_Geometry_Geometry_service
	  larsim_Simulation
	  lardataobj_Simulation
	  lardata_Utilities
	  larevt_Filters
	  lardataobj_RawData
	  lardataobj_RecoBase
	  lardata_RecoObjects
	  larreco_RecoAlg
	  larpandora_LArPandoraInterface
	  nusimdata_SimulationBase
	  ${LARRECO_LIB}
	  ${LARDATA_LIB}
	  ${ART_FRAMEWORK_CORE}
	  ${ART_FRAMEWORK_PRINCIPAL}
	  ${ART_FRAMEWORK_SERVICES_REGISTRY}
	  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	  art_Persistency_Common
	  art_Persistency_Provenance
	  art_Utilities
	  ${MF_MESSAGELOGGER}
	  ${MF_UTILITIES}
	  ${FHICLCPP}
	  ${CETLIB}
	  ${ROOT_GEOM}
	  ${ROOT_XMLIO}
	  ${ROOT_GDML}
	  ${ROOT_BASIC_LIB_LIST}
	  ${LARCV_LIB}
	  ${BOOST_LIB}
        )

install_headers()
install_fhicl()
install_source()
