import static groovy.io.FileTypes.FILES

new File('.').eachFileRecurse(FILES){
    if(it.name.endsWith('.groovy')){
        println it
    }
}
