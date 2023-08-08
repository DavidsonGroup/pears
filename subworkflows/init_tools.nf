

process initialise_tools{
    if checkIfExists('') 
    else 
        cd $projectDir/modules/flexiplex && make

    if checkIfExists('')
    else
        cd $projectDir/modules/arriba && make

}
