

process initialise_tools{
    if checkIfExists('') 
    else 
        cd $projectDir/modules/flexiplex && make

    if checkIfExists('')
    else
        cd $projectDir/modules/arriba && make

    if checkIfExists('')
	cd $projectDir/modules/STAR/source && make

}
