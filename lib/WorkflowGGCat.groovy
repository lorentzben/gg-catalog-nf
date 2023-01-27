/* groovylint-disable TrailingWhitespace */
//
// This file holds several functions specific to the workflow/ampliseq.nf in the nf-core/ampliseq pipeline
// This is based one the WorkflowAmpliseq.groovy file from nf-core/ampliseq

import groovy.text.SimpleTemplateEngine

class WorkflowGGCat {  
    //
    // Check string (String s) ends with one entry of an array of strings ("String[] extn")
    //
    public static boolean checkIfFileHasExtension(String s, String[] extn) {
        return Arrays.stream(extn).anyMatch(entry -> s.endsWith(entry));
    }
}