source s00.load_environment.sh


msmc2_Linux \
	--skipAmbiguous \
	--outFilePrefix ${WORK_DIR}/output/${RUN_NAME} \
	--I 0-2,0-3,1-2,1-3 \
	--timeSegmentPattern 1*2+25*1+1*2+1*3 \
    --outFilePrefix ${WORK_DIR}/output/${RUN_NAME} \
	${WORK_DIR}/output/multihetsep.txt 
