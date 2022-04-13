function mobs_data ()
--ADD eyes falgs
--REMOVE mvcoff ?
mobs_stats={}
mobs_stats.goblin={class="goblin", subclass=1, race="goblin",  sprite = "goblin", name="aarrrgh", face=5, lv=1, nature="humanoid", size="normal", gender=1,  age=15, motion="walking", track=2, hexes=1,
            mgt=15,enu=15,dex=15,spd=15,acu=15,sns=15,int=5,spr=5,chr=5,luk=5, --stats
            hp_base_mod = 12, sp_base_mod = 0, st_base_mod = 12,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0,moral = 30, fov = 180, nightvision = 1, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","rf","lf","body","body","body","body","body"}},
            mvcoff=1, expa=100,
            perks={},
            rezfire=0, rezcold=10, rezstatic=0, rezpoison=20, rezacid=5, rezmind=0, rezspirit=0, rezdisease=20, rezlight=0, rezdarkness=0, --resistances
			immunities = {},
			revenge_type = 0, revenge_power = 0,
            spellbook=0, warbook=0, spellnames={}, battleai="melee", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0, multiattack = 1,--battle stats
            num_unarmed=0,num_dagger=0,num_sword=3,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=0,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=0,lvl_dagger=0,lvl_sword=1,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=0,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=3, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=1, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=0,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=0,num_mysticism=0,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=5, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=0,lvl_mysticism=0,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=1, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0, lvl_leadership = 0, lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            --inventory_list={{ttxid=16,q=math.random(1,inventory_ttx[16].material),w=0,e=0,r=0,h=0,v=50},{ttxid=53,q=math.random(1,inventory_ttx[53].material),w=0,e=0,r=0,h=0,v=50},{ttxid=413,q=1,w=0,e=0,r=0,v=20},{ttxid=415,q=1,w=0,e=0,r=0,v=75},{ttxid=416,q=1,w=0,e=0,r=0,h=0,v=50}},
			inventory_list={{ttxid=16,q=math.random(1,inventory_ttx[16].material),w=0,e=0,r=0,h=0,v=50},{ttxid=53,q=math.random(1,inventory_ttx[53].material),w=0,e=0,r=0,h=0,v=50}},
			equipment={rh=1,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0,horns=0,tail=0},
			flayloot={},
			loot_gold_min=0, loot_gold_max=10
                   };
mobs_stats.rogue={class="rogue", subclass=1, race="human",  sprite = "rogue", name="nonameandnick", face=5, lv=1, nature="humanoid", size="normal",  age=25, gender=1, motion="walking",  track=1, hexes=1,
            mgt=10,enu=15,dex=20,spd=20,acu=35,sns=55,int=5,spr=5,chr=5,luk=25, --stats
            hp_base_mod = 12, sp_base_mod = 0, st_base_mod = 12,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0, moral = 30, fov = 180, nightvision = 0, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","rf","lf","body","body","body","body","body"}},
            mvcoff=1, recmel=50, recrng=50, recmag=0, expa=100,
            perks={},
            rezfire=0, rezcold=10, rezstatic=0, rezpoison=20, rezacid=5, rezmind=0, rezspirit=0, rezdisease=20, rezlight=0, rezdarkness=0, --resistances
            immunities = {},
            revenge_type = 0, revenge_power = 0,
            spellbook=0, warbook=0, spellnames={}, battleai="shooter", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0, multiattack = 1,--battle stats
            num_unarmed=0,num_dagger=0,num_sword=3,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=10,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=0,lvl_dagger=0,lvl_sword=1,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=3,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=15, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=3, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=0,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=0,num_mysticism=0,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=10, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=0,lvl_mysticism=0,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=3, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0, lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            --inventory_list={{ttxid=22,q=math.random(1,inventory_ttx[22].material),w=0,e=0,r=0,h=0,v=50},{ttxid=131,q=math.random(1,inventory_ttx[131].material),w=0,e=0,r=0,h=0,v=50},{ttxid=146,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=147,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=148,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=413,q=1,w=0,e=0,r=0,v=20},{ttxid=415,q=1,w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=1,w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=1,w=0,e=0,r=0,h=0,v=50}},
			inventory_list={{ttxid=22,q=math.random(1,inventory_ttx[22].material),w=0,e=0,r=0,h=0,v=50},{ttxid=131,q=math.random(1,inventory_ttx[131].material),w=0,e=0,r=0,h=0,v=50},{ttxid=146,q=math.random(10,30),w=0,e=0,r=0,h=0,v=50},{ttxid=147,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=148,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50}},
			equipment={rh=1,lh=0,ranged=2,ammo=3, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
			flayloot={},
			loot_gold_min=0, loot_gold_max=10
                   };
mobs_stats.mage={class="mage", subclass=1, race="human",  sprite = "rogue", name="nonameandnick", face=5, lv=1, nature="humanoid", size="normal",  age=35, gender=1, motion="walking",  track=1, hexes=1,
            mgt=5,enu=10,dex=5,spd=5,acu=5,sns=15,int=55,spr=0,chr=0,luk=25, --stats
            hp_base_mod = 12, sp_base_mod = 25, st_base_mod = 12,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0, moral = 30, fov = 180, nightvision = 0, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","rf","lf","body","body","body","body","body"}},
            mvcoff=1, recmel=50, recrng=50, recmag=0, expa=100,
            perks={},
            shader={225, 200, 0, 125},
            rezfire=0, rezcold=10, rezstatic=0, rezpoison=20, rezacid=5, rezmind=0, rezspirit=0, rezdisease=20, rezlight=0, rezdarkness=0, --resistances
            immunities = {},
            revenge_type = 0, revenge_power = 0,
            spellbook=1, spellnames={"flamearrow","fireball","fireblast","coldbeam","poisonedspit","acidburst","iceball"}, sp_limit = 0,battleai="battlemage",protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0, multiattack = 1,--battle stats
            num_unarmed=0,num_dagger=0,num_sword=3,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=0,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=0,lvl_dagger=0,lvl_sword=1,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=0,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=15, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=3, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=15,num_water=15,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=3,lvl_water=3,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=15,num_mysticism=15,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=3,lvl_mysticism=3,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0, lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
			inventory_list={{ttxid=86,q=math.random(1,inventory_ttx[86].material),w=0,e=0,r=0,h=0,v=50}},
			equipment={rh=1,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
			flayloot={},
			loot_gold_min=0, loot_gold_max=10
                   };
                   
mobs_stats.ratman={class="ratman", subclass=1, race="ratman",  sprite = "rogue", name="eeep", face=5, lv=5, nature="humanoid", size="normal",  age=15, gender=1, motion="walking",  track=1, hexes=1,
            mgt=10,enu=15,dex=25,spd=25,acu=25,sns=25,int=10,spr=5,chr=5,luk=25, --stats
            hp_base_mod = 10, sp_base_mod = 0, st_base_mod = 10,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0, moral = 20, fov = 180, nightvision = 1, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","rf","lf","body","body","body","body","body"}},
            mvcoff=1, recmel=50, recrng=50, recmag=0, expa=100,
            perks={},
            rezfire=0, rezcold=0, rezstatic=0, rezpoison=50, rezacid=0, rezmind=0, rezspirit=0, rezdisease=50, rezlight=0, rezdarkness=0, --resistances
            immunities = {},
            revenge_type = 0, revenge_power = 0,
            spellbook=0, warbook=0, spellnames={}, battleai="melee", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0, multiattack = 1,--battle stats
            num_unarmed=0,num_dagger=0,num_sword=3,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=10,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=0,lvl_dagger=0,lvl_sword=1,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=3,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=15, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=3, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=0,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=0,num_mysticism=0,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=0,lvl_mysticism=0,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0, lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            --inventory_list={{ttxid=22,q=math.random(1,inventory_ttx[22].material),w=0,e=0,r=0,h=0,v=50},{ttxid=131,q=math.random(1,inventory_ttx[131].material),w=0,e=0,r=0,h=0,v=50},{ttxid=146,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=147,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=148,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=413,q=1,w=0,e=0,r=0,v=20},{ttxid=415,q=1,w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=1,w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=1,w=0,e=0,r=0,h=0,v=50}},
			inventory_list={{ttxid=11,q=math.random(1,inventory_ttx[11].material),w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=0,w=0,e=0,r=0,h=0,v=30},{ttxid=414,q=0,w=0,e=0,r=0,h=0,v=30}},
			equipment={rh=1,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
			flayloot={},
			loot_gold_min=0, loot_gold_max=10
                   };

mobs_stats.wererat={class="wererat", subclass=1, race="ratman",  sprite = "rogue", name="eeep", face=5, lv=10, nature="humanoid", size="normal",  age=15, gender=1, motion="walking",  track=1, hexes=1,
            mgt=15,enu=20,dex=30,spd=30,acu=30,sns=25,int=10,spr=5,chr=5,luk=25, --stats
            hp_base_mod = 12, sp_base_mod = 0, st_base_mod = 10,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0, moral = 30, fov = 180, nightvision = 1, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","rf","lf","body","body","body","body","body"}},
            mvcoff=1, recmel=50, recrng=50, recmag=0, expa=100,
            perks={},
            rezfire=5, rezcold=5, rezstatic=5, rezpoison=50, rezacid=5, rezmind=5, rezspirit=5, rezdisease=50, rezlight=0, rezdarkness=0, --resistances
            immunities = {},
            revenge_type = 0, revenge_power = 0,
            spellbook=0, warbook=0, spellnames={}, battleai="melee", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0, multiattack = 1,--battle stats
            num_unarmed=0,num_dagger=0,num_sword=10,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=10,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=0,lvl_dagger=0,lvl_sword=3,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=3,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=15, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=3, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=0,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=0,num_mysticism=0,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 5, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=0,lvl_mysticism=0,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0, lvl_leadership = 2,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            --inventory_list={{ttxid=22,q=math.random(1,inventory_ttx[22].material),w=0,e=0,r=0,h=0,v=50},{ttxid=131,q=math.random(1,inventory_ttx[131].material),w=0,e=0,r=0,h=0,v=50},{ttxid=146,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=147,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=148,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=413,q=1,w=0,e=0,r=0,v=20},{ttxid=415,q=1,w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=1,w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=1,w=0,e=0,r=0,h=0,v=50}},
			inventory_list={{ttxid=11,q=math.random(1,inventory_ttx[11].material),w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=0,w=0,e=0,r=0,h=0,v=30},{ttxid=414,q=0,w=0,e=0,r=0,h=0,v=30}},
			equipment={rh=1,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
			flayloot={},
			loot_gold_min=5, loot_gold_max=25
                   };
                   
mobs_stats.assasinrat={class="assasinrat", subclass=1, race="ratman",  sprite = "rogue", name="eeep", face=5, lv=25, nature="humanoid", size="normal",  age=15, gender=1, motion="walking",  track=1, hexes=1,
            mgt=25,enu=20,dex=50,spd=50,acu=50,sns=50,int=10,spr=5,chr=5,luk=90, --stats
            hp_base_mod = 12, sp_base_mod = 0, st_base_mod = 10,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0, moral = 50, fov = 180, nightvision = 1, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","rf","lf","body","body","body","body","body"}},
            mvcoff=1, recmel=50, recrng=50, recmag=0, expa=100,
            perks={},
            rezfire=20, rezcold=20, rezstatic=20, rezpoison=50, rezacid=20, rezmind=20, rezspirit=20, rezdisease=50, rezlight=0, rezdarkness=0, --resistances
            immunities = {},
            revenge_type = 0, revenge_power = 0,
            spellbook=0, warbook=0, spellnames={}, battleai="melee", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0,  multiattack = 1,--battle stats
            num_unarmed=0,num_dagger=20,num_sword=0,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=10,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=0,lvl_dagger=5,lvl_sword=0,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=3,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=20, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=5, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=0,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=0,num_mysticism=0,num_trading=0,num_diplomacy=0, num_stealth=20, num_picklocking=0, num_traps=0, num_spothidden=20, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=0,lvl_mysticism=0,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=5, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=5, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0, lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            --inventory_list={{ttxid=22,q=math.random(1,inventory_ttx[22].material),w=0,e=0,r=0,h=0,v=50},{ttxid=131,q=math.random(1,inventory_ttx[131].material),w=0,e=0,r=0,h=0,v=50},{ttxid=146,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=147,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=148,q=math.random(1,10),w=0,e=0,r=0,h=0,v=50},{ttxid=413,q=1,w=0,e=0,r=0,v=20},{ttxid=415,q=1,w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=1,w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=1,w=0,e=0,r=0,h=0,v=50}},
			inventory_list={{ttxid=93,q=math.random(1,inventory_ttx[93].material),w=0,e=0,r=0,h=0,v=50},{ttxid=90,q=math.random(1,inventory_ttx[90].material),w=0,e=0,r=0,h=0,v=50},{ttxid=418,q=0,w=0,e=0,r=0,h=0,v=30},{ttxid=252,q=math.random(1,inventory_ttx[252].material),w=0,e=0,r=0,h=0,v=30},{ttxid=414,q=0,w=0,e=0,r=0,h=0,v=30}},
			equipment={rh=1,lh=2,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
			flayloot={},
			loot_gold_min=10, loot_gold_max=50
                   };
                   
mobs_stats.golem={class="golem", subclass=1, race="golem",  sprite = "golem", name="гркх", face=9,  age=1, lv=1, nature="golem", size="normal", gender=1, motion="walking",  track=9, hexes=1,--personalities
            mgt=20,enu=35,dex=15,spd=10,acu=10,sns=10,int=1,spr=0,chr=0,luk=0, --stats
            hp_base_mod = 25, sp_base_mod = 0, st_base_mod = 15,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0,moral = 0,fov = 180, nightvision = 0, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","rf","lf","body","body","body","body","body"}},
            mvcoff=1,  expa=250,
            perks={},
            rezfire=70, rezcold=70, rezstatic=70, rezpoison=20, rezacid=10, rezmind=100, rezspirit=100, rezdisease=100, rezlight=0, rezdarkness=0, --resistances
			immunities = {"stone","freeze"},
			revenge_type = 0, revenge_power = 0,
            spellbook=0, warbook=0, spellnames={}, battleai="melee", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=10, dr=50, block=0, amel=2,bmel=6,cmel=5, atkm=3, arng=0,brng=0,crng=0, atkr=0, multiattack = 1,--battle stats
            num_unarmed=15,num_dagger=0,num_sword=0,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=0,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=3,lvl_dagger=0,lvl_sword=0,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=0,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=0, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=0, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=0,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=0,num_mysticism=0,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=0,lvl_mysticism=0,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0,lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            inventory_list={}, 
            equipment={rh=0,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
            flayloot={},
            loot_gold_min=0, loot_gold_max=0
                   };
mobs_stats.fireelemental={class="fireelemental", subclass=1, race="elemental",  sprite = "fireelemental", name="пшш", face=9,  age=50, lv=1, nature="elemental", size="normal", gender=1, motion="levitation",  track=9, hexes=1,hitzones={{"head","rh","lh","rf","lf","body","body","body","body","body"}},
            mgt=10,enu=25,dex=15,spd=15,acu=15,sns=10,int=50,spr=0,chr=0,luk=0, --stats
            hp_base = 25, sp_base = 100, st_base = 0,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0,moral = 0,fov = 180, nightvision = 0, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},
            perks={},
            shader={225, 200, 0, 125},
            mvcoff=1,  expa=250,
            rezfire=100, rezcold=0, rezstatic=100, rezpoison=100, rezacid=100, rezmind=100, rezspirit=0, rezdisease=100, rezlight=100, rezdarkness=100, --resistances
			immunities = {"stone","freeze"},
			revenge_type = "fire", revenge_power = 5,
            spellbook=1, spellnames={"fireball","fireblast"}, battleai="battlemage", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0, block=0, amel=2,bmel=6,cmel=5, atkm=3, arng=0,brng=0,crng=0, atkr=0, multiattack = 1,--battle stats
            num_unarmed=15,num_dagger=0,num_sword=0,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=0,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=3,lvl_dagger=0,lvl_sword=0,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=0,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=0, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=0, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=20,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=5,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=100,num_mysticism=100,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=5,lvl_mysticism=5,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0,lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            inventory_list={}, 
            equipment={rh=0,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
            flayloot={},
            loot_gold_min=0, loot_gold_max=0
                   };
mobs_stats.airelemental={class="airelemental", subclass=1, race="elemental",  sprite = "airelemental", name="фью", face=9, age=50, lv=1, nature="elemental", size="normal", gender=1, motion="levitation",  track=9, hexes=1,--personalities
            mgt=10,enu=25,dex=15,spd=15,acu=15,sns=10,int=50,spr=0,chr=0,luk=0, --stats
            hp_base_mod = 12, sp_base_mod = 12, st_base_mod = 12,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0,moral = 0,fov = 180, nightvision = 0, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","body","body","body","body","tail","tail"}},
            mvcoff=1,  expa=250,
            perks={},
            shader={255,255,255, 125},
            rezfire=100, rezcold=100, rezstatic=100, rezpoison=100, rezacid=100, rezmind=100, rezspirit=0, rezdisease=100, rezlight=100, rezdarkness=100, --resistances
			immunities = {"stone","freeze"},
			revenge_type = "static", revenge_power = 5,
            spellbook=0, warbook=0, spellnames={"staticcharge","windfist","bell"}, battleai="battlemage", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0, block=0, amel=2,bmel=6,cmel=5, atkm=3, arng=0,brng=0,crng=0, atkr=0, multiattack = 1,--battle stats
            num_unarmed=15,num_dagger=0,num_sword=0,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=0,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=3,lvl_dagger=0,lvl_sword=0,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=0,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=0, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=0, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=100,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=5,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=100,num_mysticism=100,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=5,lvl_mysticism=5,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0,lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            inventory_list={}, 
            equipment={rh=0,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
            flayloot={},
            loot_gold_min=0, loot_gold_max=0
                   };
mobs_stats.waterelemental={class="waterelemental", subclass=1, race="elemental",  sprite = "waterelemental", name="буль", face=9, age=50, lv=1, nature="elemental", size="normal", gender=1, motion="walking",  track=9, hexes=1,--personalities
            mgt=10,enu=25,dex=15,spd=15,acu=15,sns=10,int=50,spr=0,chr=0,luk=0, --stats
            hp_base_mod = 12, sp_base_mod = 12, st_base_mod = 12,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0,moral = 0,fov = 180, nightvision = 0, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},hitzones={{"head","rh","lh","body","body","body","body","tail","tail"}},
            mvcoff=1,  expa=250,
            perks={},
            shader={200, 200, 255, 125},
            rezfire=0, rezcold=100, rezstatic=100, rezpoison=100, rezacid=100, rezmind=100, rezspirit=0, rezdisease=100, rezlight=100, rezdarkness=100, --resistances
			immunities = {"stone","freeze"},
			revenge_type = "cold", revenge_power = 5,
            spellbook=1, spellnames={"coldbeam","poisonedspit","acidburst","iceball"}, battleai="battlemage", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=0, dr=0, block=0, amel=2,bmel=6,cmel=5, atkm=3, arng=0,brng=0,crng=0, atkr=0, multiattack = 1,--battle stats
            num_unarmed=15,num_dagger=0,num_sword=0,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=0,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=3,lvl_dagger=0,lvl_sword=0,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=0,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=0, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=0, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=20,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=5,lvl_water=5,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=100,num_mysticism=100,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=5,lvl_mysticism=5,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0,lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            inventory_list={}, 
            equipment={rh=0,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
            flayloot={},
            loot_gold_min=0, loot_gold_max=0
                   };
mobs_stats.earthelemental={class="earthelemental", subclass=1, race="elemental",  sprite = "elemental", name="гркх", face=9, age=50, lv=1, nature="elemental", size="normal", gender=1, motion="walking",  track=9, hexes=1,--personalities
            mgt=20,enu=35,dex=15,spd=10,acu=10,sns=10,int=1,spr=0,chr=0,luk=0, --stats
            hp_base_mod = 12, sp_base_mod = 12, st_base_mod = 12,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0,moral = 0,fov = 180, nightvision = 0, claws = 0, teeth = 0, tail = 0, head = 1, arms={"rh","lh"},
            mvcoff=1,  expa=250,
            perks={},
            shader={100, 100, 100, 125},
            rezfire=70, rezcold=70, rezstatic=70, rezpoison=20, rezacid=10, rezmind=100, rezspirit=100, rezdisease=100, rezlight=0, rezdarkness=0, --resistances
			immunities = {"stone","freeze"},
			revenge_type = 0, revenge_power = 0,
            spellbook=0, warbook=0, spellnames={}, battleai="melee", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=0, lh3=0, rh4=0, lh4=0, --flags
            ac=0, dt=10, dr=50, block=0, amel=2,bmel=6,cmel=5, atkm=3, arng=0,brng=0,crng=0, atkr=0, multiattack = 1,--battle stats
            num_unarmed=15,num_dagger=0,num_sword=0,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=0,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=3,lvl_dagger=0,lvl_sword=0,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=0,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=0, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=0, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=0,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=0,num_mysticism=0,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=0,lvl_mysticism=0,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0,lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            inventory_list={},
            equipment={rh=0,lh=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
            flayloot={},
            loot_gold_min=0, loot_gold_max=0
                   };
mobs_stats.naga={class="naga", subclass=1, race="naga",  sprite = "golem", name="нья", face=5, age=50, lv=1, nature="hybrid", size="normal", gender=2, motion="walking", track=7, hexes=1,--personalities
            mgt=25,enu=40,dex=25,spd=10,acu=35,sns=10,int=20,spr=0,chr=0,luk=20, --stats
            hp_base_mod = 12, sp_base_mod = 0, st_base_mod = 12,hp_coff = 1, sp_coff = 0, st_coff = 1,
            hp_regeneration = 0, sp_regeneration = 0, st_regeneration = 0,moral = 0,fov = 180, nightvision = 0, claws = 0, teeth = 0, tail = 1, head = 1, arms={"rh","lh","rh1","lh1","rh2","lh2"},hitzones={{"head","rh","lh","rh1","lh1","rh2","lh2","body","body","body","tail","tail","tail"}},
            mvcoff=1,  expa=350,
            perks={},
            rezfire=20, rezcold=30, rezstatic=20, rezpoison=50, rezacid=20, rezmind=20, rezspirit=20, rezdisease=20, rezlight=0, rezdarkness=0, --resistances
			revenge_type = 0, revenge_power = 0,
            spellbook=0, warbook=0, spellnames={}, battleai="melee", protectionmode="none", rh=1, lh=1, rh2=1, lh2=1, rh3=1, lh3=1, rh4=1, lh4=1, --flags
            ac=0, dt=10, dr=50, block=0, amel=2,bmel=6,cmel=5, atkm=3, arng=0,brng=0,crng=0, atkr=0, multiattack = 6,--battle stats
            num_unarmed=0,num_dagger=0,num_sword=20,num_axe=0,num_crushing=0,num_flagpole=0,num_staff=0,num_bow=0,num_crossbow=0,num_throwing=0,num_firearms=0,num_blaster=0,--skills weapon num
            lvl_unarmed=0,lvl_dagger=0,lvl_sword=4,lvl_axe=0,lvl_crushing=0,lvl_flagpole=0,lvl_staff=0,lvl_bow=0,lvl_crossbow=0,lvl_throwing=0,lvl_firearms=0,lvl_blaster=0,--skills weapon lv
            num_dodging=0, num_leather=0,num_chainmail=0, num_plate=0,num_shield=0,--skills armor num
            lvl_dodging=0, lvl_leather=0,lvl_chainmail=0, lvl_plate=0,lvl_shield=0,--skills armor lv
            num_fire=0,num_water=0,num_air=0,num_earth=0,num_body=0,num_spirit=0,num_mind=0,num_light=0,num_darkness=0,--skills magic num
            lvl_fire=0,lvl_water=0,lvl_air=0,lvl_earth=0,lvl_body=0,lvl_spirit=0,lvl_mind=0,lvl_light=0,lvl_darkness=0,---skills magic lv
            num_alchemy=0,num_repair=0,num_bodybuilding=0,num_armmastery=0,num_meditation=0,num_mysticism=0,num_trading=0,num_diplomacy=0, num_stealth=0, num_picklocking=0, num_traps=0, num_spothidden=0, num_monsterid=0, num_stuffid=0, num_thievery=0, num_leadership = 0, num_regeneration = 0,num_music=0,--skills other num
            lvl_alchemy=0,lvl_repair=0,lvl_bodybuilding=0,lvl_armmastery=0,lvl_meditation=0,lvl_mysticism=0,lvl_trading=0,lvl_diplomacy=0, lvl_stealth=0, lvl_picklocking=0, lvl_traps=0, lvl_spothidden=0, lvl_monsterid=0, lvl_stuffid=0, lvl_thievery=0,lvl_leadership = 0,lvl_regeneration = 0,lvl_music = 0,lvl_astrology = 0,--skills other lvl
            inventory_list={{ttxid=16,q=math.random(1,inventory_ttx[16].material),w=0,e=0,r=0,h=0,v=50},{ttxid=16,q=math.random(1,inventory_ttx[16].material),w=0,e=0,r=0,h=0,v=50},{ttxid=16,q=math.random(1,inventory_ttx[16].material),w=0,e=0,r=0,h=0,v=50},{ttxid=16,q=math.random(1,inventory_ttx[16].material),w=0,e=0,r=0,h=0,v=50},{ttxid=16,q=math.random(1,inventory_ttx[16].material),w=0,e=0,r=0,h=0,v=50},{ttxid=16,q=math.random(1,inventory_ttx[16].material),w=0,e=0,r=0,h=0,v=50}},
			equipment={rh=1,lh=0,rh1=0,lh1=0,rh2=0,lh2=0,ranged=0,ammo=0, armor=0,head=0,boots=0,gloves=0,cloak=0,belt=0,amulet=0,ring1=0,ring2=0,ring3=0,ring4=0,ring5=0,ring6=0,art,teeth=0,horns=0,tail=0},
            flayloot={},
            loot_gold_min=0, loot_gold_max=100
                   };
end;
