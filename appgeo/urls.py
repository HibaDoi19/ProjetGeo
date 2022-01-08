from django.urls import path
from . import views
urlpatterns =[
    path('direct/',views.say_hello),
    path('resultd/',views.add),
    path('inverse/', views.say_hello2),
    path('resul/', views.inverse),
    path('acceuil/', views.acceuil),
    path('essaiA/', views.essaiA ),
    path('essaiB/',views.essaiB)

]