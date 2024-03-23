//function changeClass(idElement,myClass){
//document.getElementById(idElement).setAttribute("class", myClass);
//}
//------------------------------------
// gets the element from various browsers
function getElement(element_name){
	var elem_to_return;
	if( document.getElementById ) // this is the way the standards work
		elem_to_return = document.getElementById( element_name );
	else if( document.all ) // this is the way old msie versions work
		elem_to_return = document.all[element_name];
	else if( document.layers ) // this is the way nn4 works
		elem_to_return = document.layers[element_name];
	return elem_to_return;
}
//------------------------------------
// shows a div part and hides another
function toggleLayer( whichLayer_show, whichLayer_hide, sent_from_PDB ){
	// if we entered the function from the PDB case, we also test that the MSA was not chosen already
	var elem_show1, elem_show_2, elem_hide;

	elem_show1 = getElement(whichLayer_show);
	elem_hide = getElement(whichLayer_hide);

	if(whichLayer_show = "case_pdb_yes"){ //if the user chose MSA_no before, we hide its selection
		var MSA_no = getElement("case_MSA_no_Pseq" );
		if (MSA_no.style.display=='' || MSA_no.style.display == 'block')
			MSA_no.style.display = 'none';
	}

	elem_show1.style.display = 'block';
	elem_hide.style.display = 'none';
	// to show also the MSA part
	elem_show2 = 'none';
	if (sent_from_PDB == "true"){
		elem_show2 = getElement( "case_MSA" );
		if (elem_show2 != '' && elem_show2 != 'none')
			elem_show2.style.display = 'block';
	}
}
//------------------------------------
function show_div_name(div_id){
        var show = getElement( div_id );
        show.style.display = 'block';
}
//------------------------------------
function hide_div_name(div_id){
        var show = getElement( div_id );
        show.style.display = 'none';
}

//------------------------------------
function toggle_help(help_box_name){
	var elem, vis;
	elem = getElement(help_box_name);
	vis = elem.style;
	// if the style.display value is blank we try to figure it out here  
	if(vis.display==''&&elem.offsetWidth!=undefined&&elem.offsetHeight!=undefined)
		vis.display = (elem.offsetWidth!=0&&elem.offsetHeight!=0)?'block':'none';
	vis.display = (vis.display==''||vis.display=='block')?'none':'block';
}
//------------------------------------
function toggle_Algorithm_Options(){
	var chosen_algorithm = window.document.Guidance_form.PROGRAM.value;
	var show_bootstrap = getElement("BOOTSTRAP");
	var show_order = getElement("Output_Order");
	var show_MSA_programs_GUIDANCE = getElement("MSA_PROGRAM_GUIDANCE");
	var show_MSA_programs_GUIDANCE2_HOT = getElement("MSA_PROGRAM");
	var show_guidance_warn = getElement("GUIDANCE_WARN");
	if(chosen_algorithm == 'HoT'){
		show_bootstrap.style.display = 'none';
		show_order.style.display = 'none';
        show_guidance_warn.style.display = 'none';
        show_MSA_programs_GUIDANCE2_HOT.style.display = 'block';
		show_MSA_programs_GUIDANCE.style.display='none';
	}
	if(chosen_algorithm == 'GUIDANCE'){
		show_bootstrap.style.display = 'block';
		show_order.style.display = 'block';
		show_guidance_warn.style.display = 'block';
		show_MSA_programs_GUIDANCE2_HOT.style.display = 'none';
		show_MSA_programs_GUIDANCE.style.display='block';
	}
	if(chosen_algorithm == 'GUIDANCE2'){
		show_bootstrap.style.display = 'block';
		show_order.style.display = 'block';
		show_guidance_warn.style.display = 'none';
		show_MSA_programs_GUIDANCE2_HOT.style.display = 'block';
		show_MSA_programs_GUIDANCE.style.display='none';
	}
	//	show_MSA_programs.style.display = 'block';
}
//------------------------------------
function toggle_MSA_Options(){
  //  alert ("toggle_MSA_Options");
	var chosen_algorithm = window.document.Guidance_form.MSA_Program.value;
  //  alert ("QA1:"+chosen_algorithm);
	var show_MAFFT = getElement("MAFFT_OPTION");
    var show_PRANK = getElement("PRANK_OPTION");
    if(chosen_algorithm == 'MAFFT'){
	//	alert (chosen_algorithm);
	    show_PRANK.style.display = 'none';
      	show_MAFFT.style.display = 'block';
	}
    else if(chosen_algorithm == 'PRANK'){
		show_PRANK.style.display = 'block';
        show_MAFFT.style.display = 'none';
	}
	else
	{
		show_PRANK.style.display = 'none';
		show_MAFFT.style.display = 'none';
	}
}

//------------------------------------
function toggle_batch_file(elemId){
	//var batch_file, email;
	//batch_file = getElement("batch");
	alert ("called toggle_batch_file");
	//if (elem.value != "") {
	//	show_div_name("email");
	//}
	//else {
	//	hide_div_name("email");
	//}		
}

// use a class selector if available
const copyButtonLabel = "Copy Code";

setTimeout(() => {
	let blocks = document.querySelectorAll("pre");

	blocks.forEach((block) => {
	  // only add button if browser supports Clipboard API
	  if (navigator.clipboard) {
		let button = document.createElement("button");

		button.innerText = copyButtonLabel;
		block.appendChild(button);

		button.addEventListener("click", async () => {
		  await copyCode(block, button);
		});
	  }
	});
},10)

// setTimeout(() => {
// let blocks = document.querySelectorAll("pre");
//
// 	blocks.forEach((block) => {
// 	  // only add button if browser supports Clipboard API
// 	  if (navigator.clipboard) {
// 		let button = document.createElement("button");
//
// 		button.innerText = copyButtonLabel;
// 		block.appendChild(button);
// 		block.setAttribute("tabindex", 0);
//
// 		button.addEventListener("click", async () => {
// 		  await copyCode(block, button);
// 		});
// 	  }
// 	});
// 	},10)
async function copyCode(block, button) {
  let code = block.querySelector("code");
  let text = code.innerText;

  await navigator.clipboard.writeText(text);

  // visual feedback that task is completed
  button.innerText = "Code Copied";

  setTimeout(() => {
    button.innerText = copyButtonLabel;
  }, 700);
}


//
//
// const copyButtonLabel = "Copy Code";
//
// // use a class selector if available
// let blocks = document.querySelectorAll("pre");
//
// blocks.forEach((block) => {
//   // only add button if browser supports Clipboard API
//   if (navigator.clipboard) {
//     let button = document.createElement("button");
//
//     button.innerText = copyButtonLabel;
//     block.appendChild(button);
//
//     button.addEventListener("click", async () => {
//       await copyCode(block);
//     });
//   }
// });
//
// async function copyCode(block) {
//   let code = block.querySelector("code");
//   let text = code.innerText;
//
//   await navigator.clipboard.writeText(text);
// }