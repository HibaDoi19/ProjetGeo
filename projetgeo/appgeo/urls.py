from django.urls import path
from . import views
urlpatterns =[
    path('hello/',views.say_hello),
    path('add/',views.add),
    path('inverse/', views.say_hello2),
    path('resul/', views.inverse),
    path('acceuil/', views.acceuil)
]