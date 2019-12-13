FLAGS_NOAH := -O0 -g -autodouble -noerror_limit -FR -auto -WB -traceback
module_model_constants.o: module_model_constants.F
	$(FC) -cpp $(FFLAGS) $(FLAGS_NOAH) $(HEADER_DIRS) $<
module_sf_noahutl.o: module_sf_noahutl.F
	$(FC) -cpp $(FFLAGS) $(FLAGS_NOAH) $(HEADER_DIRS) $<
module_sf_myjsfc.o: module_sf_myjsfc.F
	$(FC) -cpp $(FFLAGS) $(FLAGS_NOAH) $(HEADER_DIRS) $<
module_sf_sfclay.o: module_sf_sfclay.F
	$(FC) -cpp $(FFLAGS) $(FLAGS_NOAH) $(HEADER_DIRS) $<
module_sf_noahlsm.o: module_sf_noahlsm.F
	$(FC) -cpp $(FFLAGS) $(FLAGS_NOAH) $(HEADER_DIRS) $<
module_sf_noahmplsm.o: module_sf_noahmplsm.F
	$(FC) -cpp $(FFLAGS) $(FLAGS_NOAH) $(HEADER_DIRS) $<

FLAGS_COMM := -O0 -g -FR -auto -WB -traceback -fpe0
nrtype.o: nrtype.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
f2008funcs.o: f2008funcs.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
nr_utility.o: nr_utility.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
expIntegral.o: expIntegral.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
spline_int.o: spline_int.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
summaFileManager.o: summaFileManager.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
multiconst.o: multiconst.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
var_lookup.o: var_lookup.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
data_types.o: data_types.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
globalData.o: globalData.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
flxMapping.o: flxMapping.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
get_ixname.o: get_ixname.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
ascii_util.o: ascii_util.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
popMetadat.o: popMetadat.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
outpt_stat.o: outpt_stat.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
time_utils.o: time_utils.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
mDecisions.o: mDecisions.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
snow_utils.o: snow_utils.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
soil_utils.o: soil_utils.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
updatState.o: updatState.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<
matrixOper.o: matrixOper.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_COMM) $(HEADER_DIRS) $<

FLAGS_SUMMA := -O0 -g -FR -auto -WB -traceback -fpe0
netcdf_util.o: netcdf_util.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
def_output.o: def_output.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
modelwrite.o: modelwrite.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
read_icond.o: read_icond.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
conv_funcs.o: conv_funcs.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
sunGeomtry.o: sunGeomtry.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
convE2Temp.o: convE2Temp.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
allocspace.o: allocspace.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
checkStruc.o: checkStruc.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
childStruc.o: childStruc.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
ffile_info.o: ffile_info.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
read_attrb.o: read_attrb.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
read_pinit.o: read_pinit.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
pOverwrite.o: pOverwrite.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
read_param.o: read_param.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
paramCheck.o: paramCheck.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
check_icond.o: check_icond.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
indexState.o: indexState.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
getVectorz.o: getVectorz.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
updateVars.o: updateVars.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
var_derive.o: var_derive.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
read_force.o: read_force.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
derivforce.o: derivforce.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
snowAlbedo.o: snowAlbedo.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
canopySnow.o: canopySnow.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
tempAdjust.o: tempAdjust.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
snwCompact.o: snwCompact.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
layerMerge.o: layerMerge.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
layerDivide.o: layerDivide.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
volicePack.o: volicePack.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
qTimeDelay.o: qTimeDelay.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
vegPhenlgy.o: vegPhenlgy.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
diagn_evar.o: diagn_evar.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
stomResist.o: stomResist.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
groundwatr.o: groundwatr.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
vegSWavRad.o: vegSWavRad.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
vegNrgFlux.o: vegNrgFlux.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
ssdNrgFlux.o: ssdNrgFlux.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
vegLiqFlux.o: vegLiqFlux.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
snowLiqFlx.o: snowLiqFlx.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
soilLiqFlx.o: soilLiqFlx.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
bigAquifer.o: bigAquifer.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
computFlux.o: computFlux.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
computResid.o: computResid.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
computJacob.o: computJacob.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
eval8summa.o: eval8summa.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summaSolve.o: summaSolve.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
systemSolv.o: systemSolv.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
varSubstep.o: varSubstep.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
opSplittin.o: opSplittin.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
coupled_em.o: coupled_em.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
run_oneHRU.o: run_oneHRU.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
run_oneGRU.o: run_oneGRU.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_type.o: summa_type.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_util.o: summa_util.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_alarms.o: summa_alarms.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_globalData.o: summa_globalData.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_defineOutput.o: summa_defineOutput.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_init.o: summa_init.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_setup.o: summa_setup.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_restart.o: summa_restart.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_forcing.o: summa_forcing.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_modelRun.o: summa_modelRun.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_writeOutput.o: summa_writeOutput.f90 
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
summa_driver.o: summa_driver.f90
	$(FC) -cpp $(FFLAGS) $(FLAGS_SUMMA) $(HEADER_DIRS) $<
