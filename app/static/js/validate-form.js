//
// validate-form.js : created by Josef Sprinzak
// function validateform() is called from Submit in index.php
// added check that captcha text was filled in correctly
//  
function validateform() {

				form = document.Guidance_form

				if (form.input_type && form.input_type.value == 'msa_set') {
					if (!form.msaSet_File || form.msaSet_File.value == "") {
						confirm("Please upload a .zip or .tar.gz archive with your set of MSAs")
						return false
					}
					// the sequences come from the MSAs: drop any hidden sequence text or file
					form.FASTA_txt.value = ""
					form.usrSeq_File.value = ""
				}
				else if ((form.usrSeq_File.value == "")&&(form.FASTA_txt.value=="")){
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
				if (form.n_alt_trees) {
					var nt = parseInt(form.n_alt_trees.value, 10)
					if (isNaN(nt) || nt < 3 || nt > 100) {
						confirm("The number of alternative MSAs to infer trees from must be between 3 and 100")
						form.n_alt_trees.focus()
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
                // reCAPTCHA is only present after cookie consent; without consent the server rate-limits instead.
                var captchaMode = document.getElementById('captcha_mode');
                if (captchaMode && captchaMode.value === 'recaptcha' && window.grecaptcha &&
                    typeof grecaptcha.getResponse === 'function' && grecaptcha.getResponse() === '') {
					confirm("The captcha box 'I'm not a robot' is required to be filled in correctly")
					return false
				}
				return true	
}