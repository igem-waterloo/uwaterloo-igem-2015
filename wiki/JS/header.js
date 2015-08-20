// scripts for header

// detects scroll for top resizing and lower nav show
$(window).scroll(function() {
  if ($(document).scrollTop() > 50) {
    $('.navbar-fixed-top').addClass('shrink');
    $('.navbar-lower').removeClass('hide-lower')
  } else {
    $('.navbar-fixed-top').removeClass('shrink');
    $('.navbar-lower').addClass('hide-lower')
  }
});

// set lower navbar
$('.navbar-lower').affix({
  offset: {top: 200}
});

// fills in lower nav with inner page links
$(document).ready(function(){
$('.link').each(function(i, obj) {
    $("#inner-page-links").append('<li><a href="#'+obj.id+'">'+obj.id+'</a></li>');
});});

// for smooth scrolling
$(function() {
  $('a[href*=#]:not([href=#])').click(function() {
    if (location.pathname.replace(/^\//,'') == this.pathname.replace(/^\//,'') && location.hostname == this.hostname) {
      var target = $(this.hash);
      target = target.length ? target : $('[name=' + this.hash.slice(1) +']');
      if (target.length) {
        $('html,body').animate({
          scrollTop: target.offset().top
        }, 500);
        return false;
      }
    }
  });
});