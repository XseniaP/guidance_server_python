//
// validate-form.js : created by Josef Sprinzak
// function validateform() is called from Submit in index.php
// added check that captcha text was filled in correctly
//  
function validateform() {

				form = document.Guidance_form
                
                
				if ((form.usrSeq_File.value == "")&&(form.FASTA_txt.value=="")){
					confirm("Please provide sequences as Text or File")
					form.usrSeq_File.focus()
					return false
				}
				if ((form.usrSeq_File.value != "")&&(form.FASTA_txt.value!="")){
					confirm("Please provide sequences as Text or File (not both)")
					form.usrSeq_File.focus()
					return false
				}
                if (form.FASTA_txt.value!=""){
                    if (form.FASTA_txt.value.length > 4500) { 
                        confirm("Pasted sequence is too large, please upload as a file")
                        form.FASTA_txt.focus()
                        return false
                    }
				}
				Seq_Type="";
				for (i=0;i<form.Seq_Type.length;i++) {
					if (form.Seq_Type[i].checked) {
						Seq_Type= form.Seq_Type[i].value;
					}
				}
				if (Seq_Type==""){
					confirm("Please indicate the sequences type: Amino Acids/Nucleotides/Codons")
					return false
				}	
				if (form.Bootstraps.value == ""){
					confirm("number of bootstrap repeats was set to default value (100)")
					form.Bootstraps.focus()
					form.Bootstraps.value=100
					return false
			    }
				if (form.Bootstraps.value > 100){
					confirm("The maximal number of bootstrap repeats is 100")
					form.Bootstraps.focus()
					form.Bootstraps.value=100
					return false
			    }
				if (form.Bootstraps.value < 2){
					confirm("The minimal number of bootstrap repeats is 2")
					form.Bootstraps.focus()
					form.Bootstraps.value=100
					return false
			    }
                if (grecaptcha.getResponse() === '') {
					confirm("The captcha box 'I'm not a robot' is required to be filled in correctly")
					return false
				}
				return true	
}