//global menu

$(function(){
      $('#global_nav, .header .btn_menu')
	  .click(function(e){
		     $('#global_nav').slideToggle();
		 });
      
      $('.btn_pop')
	  .click(function(e){
		     $('#pop_menu').slideToggle('normal');
		     $('.btn_pop').toggleClass('active');
		 });
  });